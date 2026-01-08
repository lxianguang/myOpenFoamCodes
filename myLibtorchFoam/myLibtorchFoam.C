#include "torch/torch.h"
#include "Time.H"
#include "fvCFD.H"
#include "vector.H"
#include "argList.H"
#include "IFstream.H"
#include "OFstream.H"
#include "polyPatchID.H"

#include "myFNN.H"
#include "myFunctions.H"

using namespace Foam;

int main(int argc, char *argv[]) 
{
    argList::noParallel();
    argList::validArgs.append("option");
    argList::addOption
        (
            "option",
            "word",
            "input 'training' or 'predicting'!"
        );
    Foam::argList args(argc, argv);
    word runingOption = args[1]; 

    //#include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // 设置训练随机数种子
    torch::manual_seed(36);

    Info << "================ Loading training parameters ===============" << endl;
    Info << "Loading training parameters from ./" << runTime.system() + "/trainingControlDict" << endl;
    IOdictionary trainingControlDict = readTransportProperties(mesh, runTime.system(), "trainingControlDict");  // 创建自定义字典对象
    const List<label> layerSizes   = trainingControlDict.get<List<label>>("layerSizes");     // 网络层数及每层神经元数量
    const bool   isContinue        = trainingControlDict.get<bool>("isContinue");            // 是否继续上次训练
    const size_t num_epochs        = trainingControlDict.get<scalar>("numEpochs");           // 训练次数
    const float  learningRate      = trainingControlDict.get<scalar>("learningRate");        // 学习率
    const float  schedulerRate     = trainingControlDict.get<scalar>("schedulerRate");       // 学习率衰减率
    const size_t schedulerStep     = trainingControlDict.get<scalar>("schedulerStep");       // 学习率衰减步长
    const size_t infoInterval      = trainingControlDict.get<scalar>("outputInterval");      // 训练信息输出间隔
    const vector trainingTimeList  = trainingControlDict.get<vector>("trainingTimeList");    // 初始时间、时间间隔、时刻数量
    const vector predicingTimeList = trainingControlDict.get<vector>("predictingTimeList");
    const vector maxCoords         = trainingControlDict.get<vector>("maxCoords");           // 定义流场取样范围
    const vector minCoords         = trainingControlDict.get<vector>("minCoords");
    const size_t dataInterval      = trainingControlDict.get<scalar>("dataInterval");        // 数据采样间隔

    // 定义计算时间范围和受力文件路径
    const fileName forcesDir             = "./postProcessing/forcesWing/0/force.dat";
    const fileName outputDir             = "./postProcessing/trainedFNN/";
    // const instantList timeDirs = timeSelector::select0(runTime, args);
    const instantList trainingTimeDirs   = creatTimeDirList(trainingTimeList);
    const instantList predictingTimeDirs = creatTimeDirList(predicingTimeList);
    
    // 输出训练参数
    if (!isDir(outputDir)) mkDir(outputDir);
    writingTrainingParameters(
        outputDir + "trainingControlDict",
        layerSizes,
        maxCoords,
        minCoords,
        dataInterval,
        trainingTimeList,
        predicingTimeList,
        num_epochs,
        learningRate,
        schedulerRate,
        schedulerStep,
        infoInterval
    );
    Info << "================= Creating Neural Networks =================" << endl;
    DynamicFNN myFNN(layerSizes);
    torch::nn::MSELoss criterion = torch::nn::MSELoss();  // 定义损失函数
    torch::optim::Adam optimizer(myFNN.parameters(), torch::optim::AdamOptions(learningRate));  // 定义优化器
    auto scheduler = torch::optim::StepLR(optimizer, schedulerStep, schedulerRate);  // 随着训练进行，学习率逐步减小
    Info << "FNN layerSizes      : " << layerSizes << endl;
    
    // 判断是训练还是预测
    torch::Tensor force_max, force_min;
    if (runingOption == "training" || runingOption == "t")
    {
        Info << "=============== Loading physical properties ================" << endl;
        const IOdictionary transportProperties = readTransportProperties(mesh, runTime.constant(), "transportProperties");
        const dimensionedScalar nu  = readTransportNu(transportProperties);
        const dimensionedScalar rho = readTransportRho(transportProperties);
        if(isContinue)
        {
            Info << "=============== Loading training parameters ================" << endl;
            loadingTrainedModel(
                outputDir + "trainedFNN.pt",
                myFNN
            );
        }
        Info << "================== Loading training data ===================" << endl;
        std::vector<torch::Tensor> trainingInputTensorList;
        torch::Tensor trainingOutputTensor, trainingNormalizedTensor;
        creatingTrainingTensor(
            mesh, 
            runTime,
            trainingTimeDirs,
            minCoords, 
            maxCoords, 
            dataInterval,
            forcesDir, 
            trainingInputTensorList, 
            trainingOutputTensor
        );
        normalizedTrainingTensor(
            trainingOutputTensor,
            trainingNormalizedTensor,
            force_max,
            force_min
        );
        writingPointsInformation(
            outputDir + "trainingPoints.dat",
            trainingTimeDirs,
            trainingOutputTensor
        );
        Info << "================= Training Neural Network ==================" << endl;
        myFNN_training(
            myFNN,
            num_epochs,
            infoInterval,
            outputDir + "trainingLog.dat",
            trainingInputTensorList,
            trainingNormalizedTensor,
            optimizer,
            scheduler,
            criterion
        );
        Info << "=============== Saving training parameters ================" << endl;
        savingTrainedModel(
            outputDir + "trainedFNN.pt",
            myFNN,
            outputDir + "normalized.pt",
            force_max,
            force_min
        );
    }
    else if (runingOption == "predicting" || runingOption == "p")
    {
        Info << "=============== Loading training parameters ===============" << endl;
        loadingTrainedModel(
            outputDir + "trainedFNN.pt",
            myFNN,
            outputDir + "normalized.pt",
            force_max,
            force_min
        );
        Info << "================ predicting Neural Network =================" << endl;
        std::vector<torch::Tensor> predictingInputTensorList;
        torch::Tensor predictingOutputTensor, predicted_forces;
        creatingTrainingTensor(
            mesh, 
            runTime,
            predictingTimeDirs,
            minCoords, 
            maxCoords, 
            dataInterval,
            forcesDir, 
            predictingInputTensorList, 
            predictingOutputTensor
        );
        myFNN_predicting(
            myFNN,
            predictingTimeDirs,
            predictingInputTensorList,
            predictingOutputTensor,
            force_max,
            force_min,
            predicted_forces
        );
        writingPointsInformation(
            outputDir + "predictingPoints.dat",
            predictingTimeDirs,
            predicted_forces
        );
    }
    else
    {
        FatalError << "Invalid running option! Use 'training' or 'predicting'."
                   << abort(FatalError);
    }

    return 0;
}
