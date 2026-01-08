#ifndef MYFNN_H
#define MYFNN_H

#include "torch/torch.h"

using namespace Foam;

// 定义模块实现类
struct DynamicFNN : torch::nn::Module 
{
    torch::nn::Sequential net_;

    DynamicFNN() = default;
    
    // 构造函数：适配OpenFOAM List<size_t>类型：神经元数量列表 [input_size, hidden1, hidden2, ..., output_size]
    DynamicFNN(const Foam::List<Foam::label>& layer_sizes)
    {
        TORCH_CHECK(layer_sizes.size() >= 2, "至少需要输入和输出层");
        net_ = register_module("net", torch::nn::Sequential());

        // 构建隐藏层
        for (Foam::label i = 0; i < layer_sizes.size() - 1; ++i) 
        {
            // 添加线性层，将Foam::label转换为int64_t
            net_->push_back(torch::nn::Linear(
                static_cast<int64_t>(layer_sizes[i]), 
                static_cast<int64_t>(layer_sizes[i + 1])
            ));
            // 除了最后一层外都添加ReLU激活函数
            if (i < layer_sizes.size() - 2) 
            {
                net_->push_back(torch::nn::ReLU());
            }
        }
    }

    // 归一化训练参数到[-1, 1]
    torch::Tensor normalizedCoordinates(torch::Tensor inputTensor) 
    {
        auto max_paras = inputTensor.max().item<float>();
        auto min_paras = inputTensor.min().item<float>();
        auto normalized_paras = 2.0 * (inputTensor - min_paras) / (max_paras - min_paras) - 1.0;
        return normalized_paras;
    }

    torch::Tensor forward(torch::Tensor inputTensor)
    {
        return net_->forward(normalizedCoordinates(inputTensor));
    }
};

/**
 * @brief 训练神经网络
 * @param myFNN 神经网络对象
 * @param num_epochs 训练次数
 * @param interval 信息输出间隔
 * @param trainingInputTensorList 输入数据集
 * @param trainingNormalizedTensor 输出数据集
 * @param optimizer 神经网络优化器
 * @param scheduler 动态学习率调整器
 * @param criterion 损失函数
 */
bool myFNN_training(
    DynamicFNN& myFNN,
    const size_t num_epochs,
    const size_t interval,
    const fileName& outputPath,
    std::vector<torch::Tensor>& trainingInputTensorList,
    torch::Tensor& trainingNormalizedTensor,
    torch::optim::Adam& optimizer,
    torch::optim::StepLR& scheduler,
    torch::nn::MSELoss& criterion
)
{
    autoPtr<OFstream> outputFilePtr;
    outputFilePtr.reset(new OFstream(outputPath));
    outputFilePtr() << "Variables = epoch, epoch_loss, LR" << "\n" << endl;;
    for (size_t epoch = 1; epoch < num_epochs + 1; ++epoch)
    {
        float epoch_loss = 0.0;
        myFNN.train();
        
        for (size_t timeCount = 0; timeCount < trainingInputTensorList.size(); timeCount++)
        {
            torch::Tensor train_input   = trainingInputTensorList[timeCount];
            torch::Tensor train_output  = trainingNormalizedTensor[timeCount]; 
            
            torch::Tensor output_tensor = myFNN.forward(train_input);
            torch::Tensor pred_output   = output_tensor.mean(0);
            torch::Tensor training_loss = criterion(pred_output, train_output);
            
            optimizer.zero_grad();  
            training_loss.backward();     // 反向传播
            optimizer.step();             // 优化器更新参数
            
            epoch_loss += training_loss.item<float>();
        }
    
        scheduler.step(); // 学习率调度

        // 输出训练信息
        float average_loss = epoch_loss / trainingInputTensorList.size(); // 计算平均损失
        if (epoch % interval == 0)
        {
            std::cout << "Epoch [" << epoch << "/" << num_epochs << "]" << " | "
                      << "Loss: "  << average_loss
                      << " | LR: " << optimizer.param_groups()[0].options().get_lr() << std::endl;
            outputFilePtr() << epoch        << "    " 
                            << average_loss << "    " 
                            << optimizer.param_groups()[0].options().get_lr() << endl;
        }
    }

    return true;
}

/**
 * @brief 训练神经网络
 * @param myFNN 神经网络对象
 * @param predictingInputTensorList 验证输入数据集
 * @param predictingOutputTensor 验证输出数据集
 * @param force_max 反归一化参数
 * @param force_min 反归一化参数
 * @param normalized_forces 反归一化后的预测结果
 */
bool myFNN_predicting(
    DynamicFNN& myFNN,
    const instantList& timeDirs,
    std::vector<torch::Tensor>& predictingInputTensorList,
    torch::Tensor& predictingOutputTensor,
    torch::Tensor& force_max,
    torch::Tensor& force_min,
    torch::Tensor& predicted_forces
)
{
    myFNN.eval();
    std::vector<float> outputData;
    for (size_t timeCount = 0; timeCount < predictingInputTensorList.size(); ++timeCount)
    {
        torch::Tensor validate_input    = predictingInputTensorList[timeCount];
        torch::Tensor validate_forces   = predictingOutputTensor[timeCount];
        
        torch::Tensor predicted_output1 = myFNN.forward(validate_input);
        torch::Tensor predicted_output2 = predicted_output1.mean(0);
        torch::Tensor normalized_forces = 0.5 * (predicted_output2 + 1.0) * (force_max - force_min) + force_min;
        
        outputData.push_back(normalized_forces[0].item<float>());
        outputData.push_back(normalized_forces[1].item<float>());

        std::cout << "Time: " << timeDirs[timeCount].name() 
                  << "------------" << std::endl;
        std::cout << "Predicted: (" 
                  << normalized_forces[0].item<float>() << ", " 
                  << normalized_forces[1].item<float>() << ") " << std::endl;
        std::cout << "Actual   : (" 
                  << validate_forces[0].item<float>()   << ", " 
                  << validate_forces[1].item<float>()   << ")" << std::endl;
    }
    predicted_forces = torch::from_blob(
                           outputData.data(),
                           {static_cast<int64_t>(outputData.size())},
                           torch::kFloat32
                       ).clone();
    predicted_forces = predicted_forces.reshape({-1, 2});                    

    return true;
}

/**
 * @brief 训练神经网络
 * @param outputPath1 模型文件保存路径
 * @param myFNN 神经网络对象
 * @param outputPath1 归一化参数文件保存路径
 * @param force_max 归一化参数
 * @param force_max 归一化参数
 */
bool savingTrainedModel(
    const fileName& outputPath1,
    DynamicFNN& myFNN,
    const fileName& outputPath2,
    torch::Tensor& force_max,
    torch::Tensor& force_min
)
{
    torch::serialize::OutputArchive archive;
    myFNN.save(archive);
    archive.save_to(outputPath1);
    Info << "Writing trained model to " << outputPath1 << endl;

    torch::serialize::OutputArchive normalized;
    normalized.write("force_max", force_max);
    normalized.write("force_min", force_min);
    normalized.save_to(outputPath2);
    Info << "Writing ormalization parameters to " << outputPath2 << endl;

    return true;
}

/**
 * @brief 训练神经网络
 * @param outputPath1 模型文件保存路径
 * @param myFNN 神经网络对象
 */

bool loadingTrainedModel(
    const fileName& inputPath1,
    DynamicFNN& myFNN
)
{
    if(!isFile(inputPath1))
    {
        FatalErrorIn("loadingTrainedModel")
        << "Model file " << inputPath1 << " not found ! " 
        << exit(FatalError);
    }
    torch::serialize::InputArchive archive;
    archive.load_from(inputPath1);
    myFNN.load(archive);
    Info << "Loading trained model from " << inputPath1 << endl;

    return true;
}

/**
 * @brief 训练神经网络
 * @param outputPath1 模型文件保存路径
 * @param myFNN 神经网络对象
 * @param outputPath1 归一化参数文件保存路径
 * @param force_max 归一化参数
 * @param force_max 归一化参数
 */

bool loadingTrainedModel(
    const fileName& inputPath1,
    DynamicFNN& myFNN,
    const fileName& inputPath2,
    torch::Tensor& force_max,
    torch::Tensor& force_min
)
{
    if(!isFile(inputPath1))
    {
        FatalErrorIn("loadingTrainedModel")
        << "Model file " << inputPath1 << " not found ! " 
        << exit(FatalError);
    }
    if(!isFile(inputPath2))
    {
        FatalErrorIn("loadingTrainedModel")
        << "Normalized parameter file " << inputPath2 << " not found ! " 
        << exit(FatalError);
    }
    torch::serialize::InputArchive archive;
    archive.load_from(inputPath1);
    myFNN.load(archive);
    Info << "Loading trained model from " << inputPath1 << endl;

    torch::serialize::InputArchive normalized;
    normalized.load_from(inputPath2);
    normalized.read("force_max", force_max);
    normalized.read("force_min", force_min);
    Info << "Loading Normalized parameters from " << inputPath2 << endl;

    return true;
}

#endif
