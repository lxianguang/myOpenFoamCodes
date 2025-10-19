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
    the weighted pressure source force decomposition (WPS) for a wall boundary (two wall boundaries in the system)

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
    argList::validArgs.append("boundaryName1");
    argList::validArgs.append("boundaryName2");
    argList::validArgs.append("dimensionLabel");

    // 准备选项
    argList::addOption // string variable
        (
            "boundaryName1",
            "word",
            "provide the name of boundary for the wps theory decomposition"
        );

    argList::addOption // string variable
        (
            "boundaryName2",
            "word",
            "provide the name of the other boundary"
        );
    argList::addOption // string variable
        (
            "dimensionLabel",
            "word",
            "the direction of the resultant force decomposition (x, y, z)"
        );    

    // 创建参数列表，通常已经在在 setRootCase.H 定义，所以要注释
    Foam::argList args(argc, argv);
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

    // 读取参数
    const word boundaryName1  = args[1];    // 读取力分解固壁边界名称
    const word boundaryName2  = args[2];    // 读取另一个固壁边界名称
    const word dimensionLabel = args[3];    // 读取合力分解的方向

    // 判断方向
    label forceIndex = 0;
    word phiFileName;
    if (dimensionLabel.compare("x") == 0)
	{   
        phiFileName = "Tx";
		forceIndex  = 0;
	}
	else if (dimensionLabel.compare("y") == 0)
	{
        phiFileName = "Ty";
		forceIndex  = 1;
	}
    else if (dimensionLabel.compare("z") == 0)
	{
        phiFileName = "Tz";
		forceIndex  = 2;
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
    Info << "loading density value: " << rho << endl;

    // 创建输出文件
    fileName outputDir = mesh.time().path()/"postProcessing/forceDecomposition";
    mkDir(outputDir);
    autoPtr<OFstream> outputFilePtr;
    outputFilePtr.reset(new OFstream(outputDir/"wpsDecomposition2_" + dimensionLabel + ".dat"));
    outputFilePtr() << "Variables = time, force_t, force_Q, force_f, force_a1, force_a2, force_vp1, force_vp2" << "\n" << endl;;

    // 读取不同时刻文件下的场量进行WPS受力分解
    instantList timeDirs = timeSelector::select0(runTime, args);
    forAll(timeDirs, timeI)
    {
        // 更新网格信息
        mesh.readUpdate();    

        // 读取文件时刻
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        volScalarField Phi(                         // 定义一个标量场，无需指定量纲，因为其量纲已经在相应的文件中指定了
            IOobject(
                phiFileName,                        // 指定名称
                runTime.timeName(),                 // 获取当前时间
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
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
                ),
            mesh
        );
        Info << "loading velocity field =============================" << endl;

        volVectorField acceleration(
            IOobject(
                "ddt(U)",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
                ),
            mesh
        );
        Info << "loading acceleration field =========================" << endl;

        // 计算涡量和Q准则
        volVectorField omega = fvc::curl(velocity);
        volTensorField gradU = fvc::grad(velocity);
        volVectorField curlO = fvc::curl(omega);
        volScalarField Qcriterion = 0.5*(sqr(tr(gradU)) - tr(((gradU) & (gradU))));

        // 体积分计算Q力
        scalar qForce = 0.0;
        forAll(mesh.C(), cellID){
            qForce = qForce + 2 * rho.value() *mesh.V()[cellID] * Phi[cellID] * Qcriterion[cellID];
        }
        Info << "Qcriterion vortex force value: " << qForce << endl;

        // 获取固体物面信息
        label patchID1 = mesh.boundaryMesh().findPatchID(boundaryName1);
        label patchID2 = mesh.boundaryMesh().findPatchID(boundaryName2);
        if (patchID1 == -1)
        {
            Info << "Boundary 1 '" << boundaryName1 << "' is not found." << endl;
            return 1; 
        }
        if (patchID2 == -1)
        {
            Info << "Boundary 2 '" << boundaryName2 << "' is not found." << endl;
            return 1; 
        }
        label faceINum1 = mesh.boundaryMesh()[patchID1].size();
        label faceINum2 = mesh.boundaryMesh()[patchID2].size();

        // 壁面积分
        scalar friction  = 0.0;
        scalar pressure1 = 0.0;
        scalar pressure2 = 0.0;
        scalar accforce1 = 0.0;
        scalar accforce2 = 0.0;
        for (label faceID = 0; faceID < faceINum1; faceID++){
            // 壁面信息获取
            scalar surfaceArea   = mesh.magSf().boundaryField()[patchID1][faceID];               // 壁面网格面积
            scalar surfacePhi    = Phi.boundaryField()[patchID1][faceID];                        // 壁面phi值    
            vector surfaceNormal = mesh.Sf().boundaryField()[patchID1][faceID] / surfaceArea;    // 壁面法向量
            vector surfaceOmega  = omega.boundaryField()[patchID1][faceID];                      // 壁面涡量
            vector surfaceCurl   = curlO.boundaryField()[patchID1][faceID];                      // 壁面涡量散度
            vector surfaceAcce   = acceleration.boundaryField()[patchID1][faceID];               // 壁面加速度场
            // 计算摩擦力
            friction  = friction  + rho.value() * nu.value() * surfaceArea * component(surfaceNormal ^ surfaceOmega, forceIndex);
            // 计算粘性压强力
            pressure1 = pressure1 + rho.value() * nu.value() * surfaceArea * surfacePhi * (surfaceNormal & surfaceCurl);
            // 计算加速度力
            accforce1 = accforce1 + rho.value() * surfaceArea * surfacePhi *(surfaceNormal & surfaceAcce);
        }
        for (label faceID = 0; faceID < faceINum2; faceID++){
            // 壁面信息获取
            scalar surfaceArea   = mesh.magSf().boundaryField()[patchID2][faceID];               // 壁面网格面积
            scalar surfacePhi    = Phi.boundaryField()[patchID2][faceID];                        // 壁面phi值    
            vector surfaceNormal = mesh.Sf().boundaryField()[patchID2][faceID] / surfaceArea;    // 壁面法向量
            vector surfaceCurl   = curlO.boundaryField()[patchID2][faceID];                      // 壁面涡量散度
            vector surfaceAcce   = acceleration.boundaryField()[patchID2][faceID];               // 壁面加速度场
            // 计算粘性压强力
            pressure2 = pressure2 + rho.value() * nu.value() * surfaceArea * surfacePhi * (surfaceNormal & surfaceCurl);
            // 计算加速度力
            accforce2 = accforce2 + rho.value() * surfaceArea * surfacePhi *(surfaceNormal & surfaceAcce);
        }
        // 计算合力
        scalar totalforce = qForce + accforce1 + accforce2 + friction + pressure1 + pressure2;
        Info << "Viscious friction force value: " << friction << endl;
        Info << "Viscious pressure force value: " << pressure1 + pressure2 << endl;
        Info << "Body acceleration force value: " << accforce1 + accforce2 << endl;

        // 输出数据
        outputFilePtr() << runTime.timeName() << " ";
        outputFilePtr() << totalforce << " ";
        outputFilePtr() << qForce     << " ";
        outputFilePtr() << friction   << " ";
        outputFilePtr() << accforce1  << " ";
        outputFilePtr() << accforce2  << " ";
        outputFilePtr() << pressure1  << " ";
        outputFilePtr() << pressure2  << endl;
    }
    return 0;
}
