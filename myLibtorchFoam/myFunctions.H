// myFunctions.H
#ifndef MYFUNCTIONS_H
#define MYFUNCTIONS_H

#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"
#include "OSspecific.H"

using namespace Foam;

/**
 * @brief 将训练参数写入文件（符合OpenFOAM字典格式）
 * @param outputPath 输出文件路径
 * @param layerSizes 神经网络层数及每层神经元数量
 * @param maxCoords 取样范围最大坐标
 * @param minCoords 取样范围最小坐标
 * @param trainingTimeList 训练时刻
 * @param predictingTimeList 预测时刻
 * @param num_epochs 训练次数
 * @param learningRate 学习率
 * @param schedulerRate 学习率衰减率
 * @param schedulerStep 学习率衰减步长
 * @param infoInterval 训练信息输出间隔
 */
bool writingTrainingParameters(
    const fileName& outputPath,
    const List<label>& layerSizes,
    const vector& maxCoords,
    const vector& minCoords,
    const vector& trainingTimeList,
    const vector& predictingTimeList,
    const size_t num_epochs,
    const double learningRate,
    const double schedulerRate,
    const size_t schedulerStep,
    const size_t infoInterval
)
{
    autoPtr<OFstream> outputFilePtr;
    outputFilePtr.reset(new OFstream(outputPath));
    
    // 写入OpenFOAM字典文件头
    outputFilePtr() << "/*--------------------------------*- C++ -*----------------------------------*\\" << endl;
    outputFilePtr() << "| =========                 |                                                 |" << endl;
    outputFilePtr() << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |" << endl;
    outputFilePtr() << "|  \\\\    /   O peration     | Version:  v2506                                 |" << endl;
    outputFilePtr() << "|   \\\\  /    A nd           | Website:  www.openfoam.com                      |" << endl;
    outputFilePtr() << "|    \\\\/     M anipulation  |                                                 |" << endl;
    outputFilePtr() << "\\*---------------------------------------------------------------------------*/" << endl;
    outputFilePtr() << "FoamFile" << endl;
    outputFilePtr() << "{" << endl;
    outputFilePtr() << "    version     2.0;" << endl;
    outputFilePtr() << "    format      ascii;" << endl;
    outputFilePtr() << "    class       dictionary;" << endl;
    outputFilePtr() << "    object      trainingControlDict;" << endl;
    outputFilePtr() << "}" << endl;
    outputFilePtr() << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << endl;
    outputFilePtr() << endl;
    
    // 写入神经网络层大小
    outputFilePtr() << "layerSizes         (";
    for (label i = 0; i < layerSizes.size(); ++i)
    {
        outputFilePtr() << static_cast<int64_t>(layerSizes[i]);
        if (i < layerSizes.size() - 1) outputFilePtr() << " ";
    }
    outputFilePtr() << ");" << endl;
    outputFilePtr() << endl;
    
    // 写入取样范围
    outputFilePtr() << "minCoords          (" << minCoords.x() << " " << minCoords.y() << " " << minCoords.z() << ");" << endl;
    outputFilePtr() << "maxCoords          (" << maxCoords.x() << " " << maxCoords.y() << " " << maxCoords.z() << ");" << endl;
    outputFilePtr() << endl;
    
    // 写入时间参数
    outputFilePtr() << "trainingTimeList   (" << trainingTimeList.x() << " " << trainingTimeList.y() << " " << trainingTimeList.z() << ");" << endl;
    outputFilePtr() << "predictingTimeList (" << predictingTimeList.x() << " " << predictingTimeList.y() << " " << predictingTimeList.z() << ");" << endl;
    outputFilePtr() << endl;
    
    // 写入训练参数
    outputFilePtr() << "numEpochs          " << num_epochs << ";" << endl;
    outputFilePtr() << "learningRate       " << learningRate << ";" << endl;
    outputFilePtr() << "schedulerRate      " << schedulerRate << ";" << endl;
    outputFilePtr() << "schedulerStep      " << schedulerStep << ";" << endl;
    outputFilePtr() << "outputInterval     " << infoInterval << ";" << endl;
    outputFilePtr() << endl;
    
    outputFilePtr() << "// ************************************************************************* //" << endl;
    
    Info << "layerSizes          : " << layerSizes << endl;
    Info << "minCoords           : " << minCoords << endl;
    Info << "maxCoords           : " << maxCoords << endl;
    Info << "trainingTimeList    : " << trainingTimeList << endl;
    Info << "predictingTimeList  : " << predictingTimeList << endl;
    Info << "numEpochs           : " << num_epochs << endl;
    Info << "learningRate        : " << learningRate << endl;
    Info << "schedulerRate       : " << schedulerRate << endl;
    Info << "schedulerStep       : " << schedulerStep << endl;
    Info << "outputInterval      : " << infoInterval << endl;

    Info << "Backup training parameters to " << outputPath << endl;

    return true;
}

/**
 * @brief 生成指定时间范围的时间列表
 * @param startTime 起始时间
 * @param deltaTime 时间间隔
 * @param timeNum 列表长度
 */
instantList creatTimeDirList(
    const vector& timeList
)
{
    instantList timeDirs;
    for (int i = 0; i < static_cast<int>(timeList[2]); i++) 
    {
        float time = timeList[0] + i * timeList[1];
        std::ostringstream oss;
        oss << time;  // 自动选择合适的精度
        std::string timeStr = oss.str();
        timeDirs.push_back(instant(timeStr));
    }

    return timeDirs;
}

/**
 * @brief 创建词典文件然后从中读取参数
 * @param mesh #include "createMesh.H"读取的网格
 * @param constantDir 数据文件路径
 */
IOdictionary readTransportProperties(
    const fvMesh& mesh,
    const fileName& constantDir,
    const word& fileName
)
{
    // 从transportProperties文件中读取数据
    // 先定义一个IOdictionary对象，其构造函数参数为一个IOobject对象
    IOdictionary transportProperties(
        IOobject(
            fileName,                           // 字典文件名
            constantDir,                        // 字典文件所在路径，这里为constant文件夹下
            mesh,                               // 一个objectRegistry类对象，在"createMesh.H"中定义
            IOobject::MUST_READ_IF_MODIFIED,    // 如果文件被修改，则必须重新读取
            IOobject::NO_WRITE                  // 表示文件为只读
        )
    );
    return transportProperties;
}

/**
 * @brief 文件中读取流体粘性
 * @param transportProperties 物理属性字典
 * @param nu 返回流体粘性
 */
dimensionedScalar readTransportNu(
    const IOdictionary& transportProperties
)
{
    // 从transportProperties文件中读取数据
    dimensionedScalar nu(
        "nu",                                   // 指定名称
        dimViscosity,                           // 指定scalar的量纲   
        transportProperties                     // 新版本写法，自动根据名称在字典中搜索
    );
    Info << "Loading viscosity   : " << nu.value() << endl;
    return nu;
}

/**
 * @brief 文件中读取流体粘性
 * @param transportProperties 物理属性字典
 * @param rho 返回流体密度
 */
dimensionedScalar readTransportRho(
    const IOdictionary& transportProperties
)
{
    // 从transportProperties文件中读取数据
    if (transportProperties.found("rho"))
    {
        dimensionedScalar rho(
            "rho",
            dimDensity,
            transportProperties
        );
        Info << "Loading density     : " << rho.value() << endl;
        return rho;
    }
    else
    {
        dimensionedScalar rho(
            "rho",
            dimDensity,
            1000.0
        );
        Info << "Density not found in transportProperties, using default value 1000. " << endl;
        return rho;
    }
}

/**
 * @brief 文件中读取指定时间的受力数据
 * @param forcesDir 受力数据文件路径
 * @param targetTime 目标时间
 * @param tolerance 时间匹配容差
 * @return 返回读取的受力向量
 */
vector readForcesAtTime(
    const fileName& forcesFile, 
    scalar targetTime,
    scalar tolerance = 1e-3
)
{
    // 检查文件是否存在
    if (!isFile(forcesFile))
    {
        FatalErrorIn("readForcesAtTime")
            << "Force file " << forcesFile << " not found ! " 
            << exit(FatalError);
    }
    // 读取文件
    IFstream file(forcesFile);
    if (!file.good())
    {
        FatalErrorIn("readForcesAtTime")
            << "Cannot open force file " << forcesFile << " ! "
            << exit(FatalError);
    }
    string header, line;
    file.getLine(header);  // 读取文件头
    
    //Info << "Reading force data for time: " << targetTime << endl;
    
    // 逐行读取文件
    while (file.getLine(line))
    {
        if (line.empty() || line[0] == '#') 
            continue;
            
        IStringStream iss(line);
        scalar timeValue;
        vector forceValue;
        
        iss >> timeValue >> forceValue.x() >> forceValue.y() >> forceValue.z();
        
        // 检查时间是否匹配（考虑数值误差）
        if (mag(timeValue - targetTime) < tolerance)
        {
            Info << "Read force: " << forceValue << endl;
            return forceValue;
        }
    }
    
    // 如果没有找到对应时间的力数据
    FatalErrorIn("readForcesAtTime")
        << "No force data found for time " << targetTime 
        << " in file " << forcesFile << " ! "
        << exit(FatalError);
        
    return vector::zero;
}

/**
 * @brief 读取不同时刻的流场数据和物体受力数据作为训练集
 * @param mesh 网格对象
 * @param runTime 时间对象
 * @param timeDirs 读取的时间列表
 * @param minCoords 取样范围最小坐标
 * @param maxCoords 取样范围最大坐标
 * @param forcesDir 受力数据文件路径
 * @param inputTensorList 返回的输入张量列表
 * @param outputTensor 返回的输出张量
 */
bool creatingTrainingTensor(
    fvMesh& mesh,
    Time& runTime,
    const instantList& timeDirs,
    const vector& minCoords,
    const vector& maxCoords,
    const fileName& forcesDir,
    std::vector<torch::Tensor>& inputTensorList,
    torch::Tensor& outputTensor
)
{
    // 读取不同时刻文件夹中的场值
    std::vector<float> inputData, outputData;
    forAll(timeDirs, timeI)
    {
        if (timeDirs[timeI].name() == "0") continue; // 跳过0时刻，因为0时刻没有速度和压力场
        runTime.setTime(timeDirs[timeI], timeI);     // 读取文件时刻
        mesh.readUpdate();                           // 更新网格信息
        Info << "Time = " << runTime.timeName() << endl;
        Info << "Loading pressure -------------------------------------------" << endl;
        volScalarField pressure(                     // 读取标量场，无需指定量纲，其量纲已经在相应的文件中指定
            IOobject(
                "p",                                 // 指定名称
                runTime.timeName(),                  // 获取当前时间
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh 
        );
        Info << "Loading velocity -------------------------------------------" << endl;
        volVectorField velocity(                     // 读取向量场
            IOobject(
                "U",
                runTime.timeName(), 
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
                ),
            mesh
        );  
        Info << "Loading coordinate -----------------------------------------" << endl;
        const volVectorField coordinate = mesh.C(); // 读取网格信息
        // 建立训练集
        inputData.clear();  
        forAll(mesh.cells(), cellI)
        {
            // 仅选择指定范围内的点
            if (coordinate[cellI].x() >= minCoords.x() && coordinate[cellI].x() <= maxCoords.x() &&
                coordinate[cellI].y() >= minCoords.y() && coordinate[cellI].y() <= maxCoords.y() &&
                coordinate[cellI].z() >= minCoords.z() && coordinate[cellI].z() <= maxCoords.z())
            {
                //inputData.push_back(runTime.value());
                inputData.push_back(coordinate[cellI].x());
                inputData.push_back(coordinate[cellI].y());
                inputData.push_back(velocity[cellI].x());
                inputData.push_back(velocity[cellI].y());
            }
        }
        // 创建输入 tensor
        Info << "Creating input tensor --------------------------------------" << endl;
        torch::Tensor inputTensor = torch::from_blob(
                                        inputData.data(),
                                        {static_cast<int64_t>(inputData.size())},
                                        torch::kFloat32
                                    ).clone();
        inputTensorList.push_back(inputTensor.reshape({-1, 4})); // 每一时刻的部分流场数据
        std::cout << "inputTensor  sizes  : " << inputTensorList.back().sizes()  << std::endl;
        Info << "Loading forces ---------------------------------------------" << endl;
        vector forceAtTime = readForcesAtTime(forcesDir, runTime.value(), 1e-3);  
        outputData.push_back(forceAtTime.x());
        outputData.push_back(forceAtTime.y());
        Info << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    }
    // 创建整体输出 tensor
    Info << "Creating output tensor -------------------------------------" << endl;
    outputTensor = torch::from_blob(
                        outputData.data(),
                        {static_cast<int64_t>(outputData.size())},
                        torch::kFloat32
                   ).clone();
    outputTensor = outputTensor.reshape({-1, 2});

    return true;
}

/**
 * @brief 对物体受力训练集归一化到[-1, 1]范围内
 * @param outputTensor 初始受力张量
 * @param normalizedTensor 返回的归一化后受力张量
 * @param force_max 返回的受力张量最大值
 * @param force_min 返回的受力张量最小值
 */
bool normalizedTrainingTensor(
    torch::Tensor& outputTensor,
    torch::Tensor& normalizedTensor,
    torch::Tensor& force_max,
    torch::Tensor& force_min
)
{
    Info << "Normalizinga output tensor ---------------------------------" << endl;
    force_max = std::get<0>(outputTensor.max(0));
    force_min = std::get<0>(outputTensor.min(0));
    normalizedTensor = 2.0 * (outputTensor - force_min) / (force_max - force_min) - 1.0;

    return true;
}

/**
 * @brief 将训练点的受力数据写入文件
 * @param outputPath 输出文件路径
 * @param timeDirs 时间列表
 * @param forceTensor 受力张量
 */
bool writingPointsInformation(
    const fileName& outputPath,
    const instantList& timeDirs,
    const torch::Tensor& forceTensor
)
{
    Info << "Writing training points to " << outputPath << endl;
    autoPtr<OFstream> outputFilePtr;
    outputFilePtr.reset(new OFstream(outputPath));
    outputFilePtr() << "Variables = t, f_x, f_y" << "\n" << endl;;
    for(int timeCount = 0; timeCount < timeDirs.size(); timeCount++)
    {
        if (timeDirs[timeCount].name() == "0") continue; // 跳过0时刻
        outputFilePtr() << timeDirs[timeCount].value() << "    "
                        << forceTensor[timeCount][0].item<float>() << "    "
                        << forceTensor[timeCount][1].item<float>() << endl;
    }

    return true;
}

/**
 * @brief 生成均匀网格的blockMeshDict字典文件
 */
bool createBlockMeshDict(
    const Time& runTime,
    const vector& domainMin,
    const vector& domainMax,
    const Vector<label>& cellCount
)
{
    Info << "Uniform mesh range : " << domainMin << " to " << domainMax << endl;
    Info << "Uniform mesh cells : " << cellCount.x() << " x " << cellCount.y() << " x " << cellCount.z()
         << " = " << cellCount.x()*cellCount.y()*cellCount.z() << endl;
    Info << "Writing blockMeshDict --------------------------------------" << endl;

    // 确保系统目录存在
    runTime.system();
    if(!isDir("./uniformMesh"))
    {
        mkDir("./uniformMesh");
    }
    cp(runTime.system(), "./uniformMesh/system");
    cp(runTime.constant(), "./uniformMesh/constant");
    // 创建blockMeshDict文件路径
    fileName dictPath = "./uniformMesh/system/blockMeshDict";

    // 使用OFstream创建文件
    OFstream os(dictPath);
    
    // 写入文件头
    os << "/*--------------------------------*- C++ -*----------------------------------*\\" << endl;
    os << "| =========                   |                                                |" << endl;
    os << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox          |" << endl;
    os << "|  \\\\    /   O peration     | Version:  v2406                                |" << endl;
    os << "|   \\\\  /    A nd           | Website:  www.openfoam.com                     |" << endl;
    os << "|    \\\\/     M anipulation  |                                                |" << endl;
    os << "\\*---------------------------------------------------------------------------*/" << endl;
    os << "FoamFile" << endl;
    os << "{" << endl;
    os << "    version     2.0;" << endl;
    os << "    format      ascii;" << endl;
    os << "    class       dictionary;" << endl;
    os << "    object      blockMeshDict;" << endl;
    os << "}" << endl;
    os << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << endl;
    os << endl;
    
    // 定义转换比例
    os << "scale 1.0;" << endl;
    os << endl;
    
    // 定义顶点
    os << "vertices" << endl;
    os << "(" << endl;
    os << "    (" << domainMin.x() << " " << domainMin.y() << " " << domainMin.z() << ")  // 0" << endl;
    os << "    (" << domainMax.x() << " " << domainMin.y() << " " << domainMin.z() << ")  // 1" << endl;
    os << "    (" << domainMax.x() << " " << domainMax.y() << " " << domainMin.z() << ")  // 2" << endl;
    os << "    (" << domainMin.x() << " " << domainMax.y() << " " << domainMin.z() << ")  // 3" << endl;
    os << "    (" << domainMin.x() << " " << domainMin.y() << " " << domainMax.z() << ")  // 4" << endl;
    os << "    (" << domainMax.x() << " " << domainMin.y() << " " << domainMax.z() << ")  // 5" << endl;
    os << "    (" << domainMax.x() << " " << domainMax.y() << " " << domainMax.z() << ")  // 6" << endl;
    os << "    (" << domainMin.x() << " " << domainMax.y() << " " << domainMax.z() << ")  // 7" << endl;
    os << ");" << endl;
    os << endl;
    
    // 定义块
    os << "blocks" << endl;
    os << "(" << endl;
    os << "    hex (0 1 2 3 4 5 6 7)" << endl;
    os << "    (" << cellCount.x() << " " << cellCount.y() << " " << cellCount.z() << ")" << endl;
    os << "    simpleGrading (1 1 1)" << endl;
    os << ");" << endl;
    os << endl;
    
    // 定义边界
    os << "boundary" << endl;
    os << "(" << endl;
    os << "    inlet" << endl;
    os << "    {" << endl;
    os << "        type patch;" << endl;
    os << "        faces" << endl;
    os << "        (" << endl;
    os << "            (0 4 7 3)" << endl;
    os << "        );" << endl;
    os << "    }" << endl;
    os << "    outlet" << endl;
    os << "    {" << endl;
    os << "        type patch;" << endl;
    os << "        faces" << endl;
    os << "        (" << endl;
    os << "            (1 2 6 5)" << endl;
    os << "        );" << endl;
    os << "    }" << endl;
    os << "    walls" << endl;
    os << "    {" << endl;
    os << "        type patch;" << endl;
    os << "        faces" << endl;
    os << "        (" << endl;
    os << "            (0 1 5 4)" << endl;
    os << "            (2 3 7 6)" << endl;
    os << "            (0 3 2 1)" << endl;
    os << "            (4 5 6 7)" << endl;
    os << "        );" << endl;
    os << "    }" << endl;
    os << ");" << endl;
    os << endl;
    
    os << "// ************************************************************************* //" << endl;
    return true;
}

/**
 * @brief 使用OpenFOAM的system函数执行blockMesh命令
 */
bool executeBlockMesh()
{
    Info << "Running blockMesh ------------------------------------------" << endl;
    
    // 使用OpenFOAM的system函数
    int result = system("cd uniformMesh && blockMesh > log.blockMesh &");
    if (result != 0)
    {
        FatalErrorIn("executeBlockMesh")
            << "blockMesh failed with exit code: " << result
            << exit(FatalError);
    }
    system("cd ..");
    return true;
}

#endif
