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
    argList::validArgs.append("boundaryName");

    // 准备选项
    argList::addOption // string variable
        (
            "boundaryName",
            "char",
            "provide the name of boundary for the WPS theory decomposition"
        );

    // 创建参数列表，通常已经在在 setRootCase.H 定义，所以要注释
    Foam::argList args(argc, argv);
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

    // 读取参数
    const word boundaryName = args[1];         // 读取固壁边界名称

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
    outputFilePtr.reset(new OFstream(outputDir/"wpsDecomposition.dat"));
    outputFilePtr() << "# Time    total_force_x(V2) total_force_y(V3) total_force_z(V4)    \
    viscous_force_x(V5) viscous_force_y(V6) viscous_force_z(V7)    \
    Qcitation_force_x(V8) Qcitation_force_y(V9) Qcitation_force_z(V10)    \
    viscous_pressure_force_x(V11) viscous_pressure_force_y(V12) viscous_pressure_force_z(V13)    \
    acceleration_force_x(V14) acceleration_force_y(V15) acceleration_force_y(V16)" << "\n" << endl;

    // 读取不同时刻文件下的场量进行WPS受力分解
    instantList timeDirs = timeSelector::select0(runTime, args);
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        mesh.readUpdate();                          // Check for new mesh
        Info<< "Time = " << runTime.timeName() << endl;

        // 读取标量场，不需要查找关键字
        volScalarField Q(                           // 定义一个标量场，无需指定量纲，因为其量纲已经在相应的文件中指定了
            IOobject(
                "Q",                                // 指定名称
                runTime.timeName(),                 // 获取当前时间
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh 
        );
        Info << "loading Q citation field" << endl;

        volScalarField Phix(
            IOobject(
                "Tx",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh 
        );

        volScalarField Phiy(
            IOobject(
                "Ty",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh 
        );

        volScalarField Phiz(
            IOobject(
                "Tz",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh 
        );
        Info << "loading phi field" << endl;

        // 读取向量场
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
        Info << "loading acceleration field" << endl;

        volVectorField omega(
            IOobject(
                "vorticity",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
                ),
            mesh
        );
        Info << "loading vorticity field" << endl;
        
        // 流场体积分计算涡力
        const scalarField field_f_Q_x = 2 * rho.value() * Phix.field() * Q.field();
        const scalarField field_f_Q_y = 2 * rho.value() * Phiy.field() * Q.field();
        const scalarField field_f_Q_z = 2 * rho.value() * Phiz.field() * Q.field();
        const scalar value_f_Q_x = gSum(mesh.V() * field_f_Q_x);
        const scalar value_f_Q_y = gSum(mesh.V() * field_f_Q_y);
        const scalar value_f_Q_z = gSum(mesh.V() * field_f_Q_z);
        Info << "Vortex force in x direction value: " << value_f_Q_x << endl;
        Info << "Vortex force in y direction value: " << value_f_Q_y << endl;
        Info << "Vortex force in z direction value: " << value_f_Q_z << endl;

        // 获取固体物面信息
        polyPatchID topPatch(boundaryName, mesh.boundaryMesh());
        if (!topPatch.active())
        {
            FatalError
                << "Patch name " << boundaryName << " not found."
                << abort(FatalError);
        }

        label patchID = topPatch.index();
        const surfaceVectorField normal = - mesh.Sf()/mesh.magSf();         // 法向量场(从物面指向流体)
        const surfaceScalarField area   = mesh.magSf();                     // 网格面积场
        const vectorField surfaceNormal = normal.boundaryField()[patchID];  // 壁面法向量
        const vectorField surfaceOmega  = omega.boundaryField()[patchID];   // 壁面涡量
        const scalarField surfaceArea   = area.boundaryField()[patchID];    // 壁面网格面积
        const scalarField surfacePhix   = Phix.boundaryField()[patchID];    // 壁面Phi值
        const scalarField surfacePhiy   = Phiy.boundaryField()[patchID];
        const scalarField surfacePhiz   = Phiz.boundaryField()[patchID];
        
        // 物面积分计算摩擦力
        const vectorField field_f_V   = rho.value() * nu.value() * (surfaceOmega ^ surfaceNormal);
        const scalarField field_f_V_x = surfaceArea * field_f_V.component(0); 
        const scalarField field_f_V_y = surfaceArea * field_f_V.component(1); 
        const scalarField field_f_V_z = surfaceArea * field_f_V.component(2); 
        const scalar value_f_V_x = gSum(field_f_V_x);
        const scalar value_f_V_y = gSum(field_f_V_y);
        const scalar value_f_V_z = gSum(field_f_V_z);
        Info << "Viscious force in x direction value: " << value_f_V_x << endl;
        Info << "Viscious force in y direction value: " << value_f_V_y << endl;
        Info << "Viscious force in z direction value: " << value_f_V_z << endl;

        // 物面积分计算粘性压强力
        const volVectorField curlOmega = fvc::curl(omega);
        const vectorField surfaceCurlOmega = rho.value() * nu.value() * curlOmega.boundaryField()[patchID];
        const scalarField field_f_VP_x = - surfaceArea * surfacePhix * (surfaceNormal & surfaceCurlOmega); 
        const scalarField field_f_VP_y = - surfaceArea * surfacePhiy * (surfaceNormal & surfaceCurlOmega); 
        const scalarField field_f_VP_z = - surfaceArea * surfacePhiz * (surfaceNormal & surfaceCurlOmega); 
        const scalar value_f_VP_x = gSum(field_f_VP_x);
        const scalar value_f_VP_y = gSum(field_f_VP_y);
        const scalar value_f_VP_z = gSum(field_f_VP_z);
        Info << "Viscious pressure force in x direction value: " << value_f_VP_x << endl;
        Info << "Viscious pressure force in y direction value: " << value_f_VP_y << endl;
        Info << "Viscious pressure force in z direction value: " << value_f_VP_z << endl;

        // 物面积分计算加速度力
        const vectorField surfaceAcceleration = rho.value() * acceleration.boundaryField()[patchID];        // 壁面加速度场
        const scalarField field_f_A_x = - surfaceArea * surfacePhix * (surfaceNormal & surfaceAcceleration); 
        const scalarField field_f_A_y = - surfaceArea * surfacePhiy * (surfaceNormal & surfaceAcceleration); 
        const scalarField field_f_A_z = - surfaceArea * surfacePhiz * (surfaceNormal & surfaceAcceleration); 
        const scalar value_f_A_x = gSum(field_f_A_x);
        const scalar value_f_A_y = gSum(field_f_A_y);
        const scalar value_f_A_z = gSum(field_f_A_z);
        Info << "Acceleration force in x direction value: " << value_f_A_x << endl;
        Info << "Acceleration force in y direction value: " << value_f_A_y << endl;
        Info << "Acceleration force in z direction value: " << value_f_A_z << endl;

        // 计算合力
        const scalar force_t_x = value_f_Q_x + value_f_VP_x + value_f_A_x + value_f_V_x;
        const scalar force_t_y = value_f_Q_y + value_f_VP_y + value_f_A_y + value_f_V_y;
        const scalar force_t_z = value_f_Q_z + value_f_VP_z + value_f_A_z + value_f_V_z;

        // 输出数据
        outputFilePtr() << runTime.timeName() << "\t";
        outputFilePtr() << force_t_x    << " ";
        outputFilePtr() << force_t_y    << " ";
        outputFilePtr() << force_t_z    << "\t";
        outputFilePtr() << value_f_V_x  << " ";
        outputFilePtr() << value_f_V_y  << " ";
        outputFilePtr() << value_f_V_z  << "\t";
        outputFilePtr() << value_f_Q_x  << " ";
        outputFilePtr() << value_f_Q_y  << " ";
        outputFilePtr() << value_f_Q_z  << "\t";
        outputFilePtr() << value_f_VP_x << " ";
        outputFilePtr() << value_f_VP_y << " ";
        outputFilePtr() << value_f_VP_z << "\t";
        outputFilePtr() << value_f_A_x  << " ";
        outputFilePtr() << value_f_A_y  << " ";
        outputFilePtr() << value_f_A_z  << endl;
    }
    return 0;
}
