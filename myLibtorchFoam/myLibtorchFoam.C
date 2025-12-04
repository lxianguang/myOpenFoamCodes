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

    Info << "==================== Setting parameters ====================" << endl;
    // 定义计算时间范围和受力文件路径
    // instantList timeDirs = timeSelector::select0(runTime, args);
    // instantList trainingTimeDirs({instant("20"), instant("20.5"), instant("21"), instant("21.5"), instant("22"), instant("22.5")});
    const instantList trainingTimeDirs   = creatTimeDirList(20.0, 0.40, 11);
    const instantList validatingTimeDirs = creatTimeDirList(20.0, 0.10, 51);
    const fileName forcesDir = "./postProcessing/forcesWing/0/force.dat";
    const fileName outputDir = "./postProcessing/trainedFNN/";

    // 定义流场取样范围
    const vector minCoords(-1.0, -1.5, 0.0);
    const vector maxCoords( 3.5,  1.5, 1.0);

    // 定义神经网络结构及训练参数
    const std::vector<int64_t> layerSizes = {4, 32, 32, 32, 32, 2}; // 网络层数及每层神经元数量
    const size_t num_epochs    = 4000;    // 训练次数
    const double learningRate  = 4e-3;    // 学习率
    const double schedulerRate = 0.9;     // 学习率衰减率
    const size_t schedulerStep = 200;     // 学习率衰减步长
    const size_t infoInterval  = 20;      // 训练信息输出间隔

    // 输出参数信息
    if (!isDir(outputDir)) mkDir(outputDir);
    writingTrainingParameters(
        outputDir + "trainingParameters.dat",
        layerSizes,
        num_epochs,
        learningRate,
        schedulerRate,
        schedulerStep,
        infoInterval,
        minCoords,
        maxCoords
    );
    Info << "================= Creating Neural Networks =================" << endl;
    DynamicFNN myFNN(layerSizes);
    torch::nn::MSELoss criterion = torch::nn::MSELoss();  // 定义损失函数
    torch::optim::Adam optimizer(myFNN.parameters(), torch::optim::AdamOptions(learningRate));  // 定义优化器
    auto scheduler = torch::optim::StepLR(optimizer, schedulerStep, schedulerRate);  // 随着训练进行，学习率逐步减小
    Info << "FNN layerSizes      : " << layerSizes << endl;
    
    // 判断是训练还是预测
    if (runingOption == "training" || runingOption == "t")
    {
        Info << "=============== Loading physical properties ================" << endl;
        const IOdictionary transportProperties = readTransportProperties(mesh, runTime.constant());
        const dimensionedScalar nu  = readTransportNu(transportProperties);
        const dimensionedScalar rho = readTransportRho(transportProperties);
        Info << "================== Loading training data ===================" << endl;
        std::vector<torch::Tensor> trainingInputTensorList;
        torch::Tensor trainingOutputTensor, trainingNormalizedTensor;
        torch::Tensor force_max, force_min;
        creatingTrainingTensor(
            mesh, 
            runTime,
            trainingTimeDirs,
            minCoords, 
            maxCoords, 
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
        torch::Tensor force_max, force_min;
        loadingTrainedModel(
            outputDir + "trainedFNN.pt",
            myFNN,
            outputDir + "normalized.pt",
            force_max,
            force_min
        );
        Info << "================ Validating Neural Network =================" << endl;
        std::vector<torch::Tensor> validatingInputTensorList;
        torch::Tensor validatingOutputTensor, predicted_forces;
        creatingTrainingTensor(
            mesh, 
            runTime,
            validatingTimeDirs,
            minCoords, 
            maxCoords, 
            forcesDir, 
            validatingInputTensorList, 
            validatingOutputTensor
        );
        myFNN_validating(
            myFNN,
            validatingTimeDirs,
            validatingInputTensorList,
            validatingOutputTensor,
            force_max,
            force_min,
            predicted_forces
        );
        writingPointsInformation(
            outputDir + "validatingPoints.dat",
            validatingTimeDirs,
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
