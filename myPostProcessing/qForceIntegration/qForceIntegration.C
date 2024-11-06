/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 by xianGuang Luo.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    force decompose

Group
    postProcessing

Description
    the weighted pressure source force decomposition (WPS) for a single wall boundary (single connected computational domain)
    this code is used to verify the convergence of the integration region

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "fvCFD.H"
#include "argList.H"
#include "OFstream.H"
#include "polyPatchID.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // 准备参数列表
    argList::noParallel();
    argList::validArgs.append("dimensionLabel");
    argList::validArgs.append("timeStep");
    argList::validArgs.append("deltaInterval");

    // 准备选项
    argList::addOption // string variable
        (
            "dimensionLabel",
            "word",
            "the direction of the resultant force decomposition (x, y, z)"
        );   

    argList::addOption // scalar variable
        (
            "timeStep",
            "scalar",
            "provide the time for the WPS theory decomposition");
    argList::addOption // string variable
        (
            "deltaInterval",
            "scalar",
            "provide the interval for vortex force integration"
        );

    // 创建参数列表，通常已经在在 setRootCase.H 定义，所以要注释
    Foam::argList args(argc, argv);
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

    // 读取参数
    const word dimensionLabel = args[1];                 // 读取合力分解的方向
    const scalar timeStep = args.get<scalar>(2);        // 读取合力分解时刻
    const scalar deltaInterval = args.get<scalar>(3);   // 读取沿流向积分间隔

    // 判断方向
    word phiFileName;
    if (dimensionLabel.compare("x") == 0)
	{   
        phiFileName = "Tx";
	}
	else if (dimensionLabel.compare("y") == 0)
	{
        phiFileName = "Ty";
	}
    else if (dimensionLabel.compare("z") == 0)
	{
        phiFileName = "Tz";
	}
    else
	{
		FatalError
                << "Dimension input " << dimensionLabel << " is illegal."
                << abort(FatalError);
	}

    //#include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"    

    // 从transportProperties文件中读取数据
    // 先定义一个IOdictionary对象，其构造函数参数为一个IOobject对象
    IOdictionary transportProperties(
        IOobject(
            "transportProperties",              // 字典文件名
            runTime.constant(),                 // 字典文件所在路径，这里为constant文件夹下
            mesh,                               // 一个objectRegistry类对象，这里没什么用
            IOobject::MUST_READ_IF_MODIFIED,    // 如果文件被修改，则必须重新读取
            IOobject::NO_WRITE                  // 表示文件为只读
        )
    );
    Info << "loading data in transportProperties" << endl;

    // 定义一个dimensionedScalar变量nu
    // nu有量纲，其量纲dimViscosity等同于(0,2,-1,0,0,0,0),单位为m2/s
    dimensionedScalar nu(
        "nu",                                   // 指定名称
        dimViscosity,                           // 指定scalar的量纲   
        transportProperties                     // com新版本写法，自动根据名称在字典中搜索
    );
    Info << "loading viscosity value: " << nu << endl;

    dimensionedScalar rho(
        "rho",
        dimDensity,
        transportProperties
    );
    Info << "loading   density value: " << rho << endl;

    // 读取指定时刻文件下的场量进行WPS受力分解
    mesh.readUpdate();                          // 更新网格
    Info<< "force decomposition time : " << timeStep << endl;

    // 读取标量场，不需要查找关键字
    volScalarField Phi(                         // 定义一个标量场，无需指定量纲，因为其量纲已经在相应的文件中指定了
        IOobject(
            phiFileName,                        // 指定名称
            runTime.timeName(timeStep),         // 获取当前时间
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh 
    );
    Info << "loading phi field ==================================" << endl;

    // 读取向量场
    volVectorField velocity(
        IOobject(
            "U",
            runTime.timeName(timeStep),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
            ),
        mesh
    );
    Info << "loading velocity field =============================" << endl;

    // 计算涡量和Q准则
    //const volVectorField omega = fvc::curl(velocity);
    const volTensorField gradU = fvc::grad(velocity);
    const volScalarField Q     = 0.5*(sqr(tr(gradU)) - tr(((gradU) & (gradU))));

    // 获取网格横向坐标范围
    const scalar xLableMin = min(mesh.points().component(0));
    const scalar xLableMax = max(mesh.points().component(0));

    // 创建输出文件
    fileName outputDir = mesh.time().path()/"postProcessing/forceDecomposition";
    mkDir(outputDir);
    autoPtr<OFstream> outputFilePtr;
    outputFilePtr.reset(new OFstream(outputDir/"wpsConvergence_" + dimensionLabel + "_" + std::to_string(timeStep).substr(0, 5) + ".dat"));
    outputFilePtr() << "# Force decomposition at time : " << timeStep << endl;
    outputFilePtr() << "# The x coordinate range of the integration region is : ";
    outputFilePtr() << "(" << xLableMin << " : " << deltaInterval << " : " << xLableMax << ")" << "\n" << endl;
    outputFilePtr() << "Variables = x, total_q_force, positive_q_force, negative_q_force" << "\n" << endl;

    // 流场体积分计算涡力
    const scalarField field_f_Q_p = 2 * rho.value() * Phi.field() * 0.5 * (Q.field() + mag(Q.field()));
    const scalarField field_f_Q_n = 2 * rho.value() * Phi.field() * 0.5 * (Q.field() - mag(Q.field()));
    //const scalar value_f_Q_p = gSum(mesh.V() * field_f_Q_p);

    for( scalar xlable = xLableMin; xlable <= xLableMax; xlable = xlable + deltaInterval )
    {
        scalar value_f_Q_p = 0.0;
        scalar value_f_Q_n = 0.0;
        // 控制积分范围
        for (label cellI = 0; cellI < mesh.C().size(); cellI++){
            if (mesh.C()[cellI].component(0) <= xlable + deltaInterval){
                value_f_Q_p = value_f_Q_p + mesh.V()[cellI] * field_f_Q_p[cellI];
                value_f_Q_n = value_f_Q_n + mesh.V()[cellI] * field_f_Q_n[cellI];
            }
        }

        scalar value_f_Q_t = value_f_Q_p + value_f_Q_n;
        
        Info << "the integral coordinate: " << xlable << endl;
        //Info << "Q-criterion force value: " << value_f_Q_t << endl;

        // 输出数据
        outputFilePtr() << xlable       << "\t";
        outputFilePtr() << value_f_Q_t  << " ";
        outputFilePtr() << value_f_Q_p  << " ";
        outputFilePtr() << value_f_Q_n  << endl;
    }
    return 0;
}
