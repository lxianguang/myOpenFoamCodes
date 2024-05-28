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
    Obtain the spanwise distribution of boundary forces

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
    argList::validArgs.append("gridsInSpanwise");

    // 准备选项
    argList::addOption // string variable
        (
            "boundaryName",
            "char",
            "provide the name of boundary for the force calculate"
        );

    argList::addOption // scalar variable
        (
            "gridsInSpanwise",
            "int",
            "Number of uniform grids in the spanwise direction");

    // 创建参数列表，通常已经在在 setRootCase.H 定义，所以要注释
    Foam::argList args(argc, argv);
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

    // 读取参数
    const word boundaryName = args[1];         // 读取固壁边界名称
    const scalar gridsInSpanwise = args.get<scalar>(2);
    
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

    // 获取网格展向坐标范围
    const scalar zLableMin = min(mesh.points().component(2));
    const scalar zLableMax = max(mesh.points().component(2));
    const scalar deltaZ    = (zLableMax - zLableMin)/gridsInSpanwise;

    // 读取文件夹内时间刻列表
    instantList timeDirs = timeSelector::select0(runTime, args);

    // 创建输出文件
    fileName outputDir = mesh.time().path()/"postProcessing/forceDecomposition";
    mkDir(outputDir);
    autoPtr<OFstream> outputFilePtr;
    outputFilePtr.reset(new OFstream(outputDir/"forceInSpanwise.dat"));
    outputFilePtr() << "# Force decomposition at boundary : " << boundaryName << endl;
    outputFilePtr() << "# The z coordinate range of the integration region is : ";
    outputFilePtr() << "(" << zLableMin << " : " << deltaZ << " : " << zLableMax << ")" << "\n" << endl;
    outputFilePtr() << "Variables = t, z, total_force_x, total_force_y, total_force_z, ";
    outputFilePtr() << "pressure_force_x, pressure_force_y, pressure_force_z, ";
    outputFilePtr() << "viscous_force_x, viscous_force_y, viscous_force_z" << endl;
    outputFilePtr() << "Zone I = " << int(gridsInSpanwise) <<", J = "<< timeDirs.size() << ", f = point" << "\n" << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        mesh.readUpdate();                          // 更新网格
        Info<< "Time = " << runTime.timeName() << endl;

        // 读取标量场，不需要查找关键字
        volScalarField p(                           // 定义一个标量场，无需指定量纲，因为其量纲已经在相应的文件中指定了
            IOobject(
                "p",                                // 指定名称
                runTime.timeName(),                 // 获取当前时间
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh 
        );
        Info << "loading pressure field" << endl;

        // 读取向量场
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
        const surfaceVectorField meshcf = mesh.Cf();                        // 网格面中心坐标
        const surfaceScalarField area   = mesh.magSf();                     // 网格面积场
        const vectorField surfaceNormal = normal.boundaryField()[patchID];  // 壁面法向量
        const vectorField surfacemeshcf = meshcf.boundaryField()[patchID];  // 壁面中心坐标
        const vectorField surfaceOmega  = omega.boundaryField()[patchID];   // 壁面涡量
        const scalarField surfacep      = p.boundaryField()[patchID];       // 壁面压力
        const scalarField surfaceArea   = area.boundaryField()[patchID];    // 壁面网格面积

        for( scalar zlable = zLableMin + 0.5 * deltaZ; zlable <= zLableMax; zlable = zlable + deltaZ )
        {
            // 物面积分计算压强力
            const vectorField field_f_P = - rho.value() * surfacep * surfaceNormal;
            const scalarField field_f_P_x = surfaceArea * field_f_P.component(0); 
            const scalarField field_f_P_y = surfaceArea * field_f_P.component(1); 
            const scalarField field_f_P_z = surfaceArea * field_f_P.component(2); 
            scalar value_f_P_x = 0.0;
            scalar value_f_P_y = 0.0;
            scalar value_f_P_z = 0.0;

            // 物面积分计算摩擦力
            const vectorField field_f_V   = rho.value() * nu.value() * (surfaceOmega ^ surfaceNormal);
            const scalarField field_f_V_x = surfaceArea * field_f_V.component(0); 
            const scalarField field_f_V_y = surfaceArea * field_f_V.component(1); 
            const scalarField field_f_V_z = surfaceArea * field_f_V.component(2); 
            scalar value_f_V_x = 0.0;
            scalar value_f_V_y = 0.0;
            scalar value_f_V_z = 0.0;

           // 控制积分范围
            for (label cellI = 0; cellI <= surfacemeshcf.size(); cellI++){
                if (fabs(surfacemeshcf[cellI].component(2) - zlable) < 0.25 * deltaZ){
                    value_f_P_x = value_f_P_x + field_f_P_x[cellI];
                    value_f_P_y = value_f_P_y + field_f_P_y[cellI];
                    value_f_P_z = value_f_P_z + field_f_P_z[cellI];
                    value_f_V_x = value_f_V_x + field_f_V_x[cellI];
                    value_f_V_y = value_f_V_y + field_f_V_y[cellI];
                    value_f_V_z = value_f_V_z + field_f_V_z[cellI];
                }
            }

            // 计算合力
            const scalar force_t_x = value_f_P_x + value_f_V_x;
            const scalar force_t_y = value_f_P_y + value_f_V_y;
            const scalar force_t_z = value_f_P_z + value_f_V_z;

            // 输出数据
            outputFilePtr() << runTime.timeName() << "\t";
            outputFilePtr() << zlable       << " ";
            outputFilePtr() << force_t_x    << " ";
            outputFilePtr() << force_t_y    << " ";
            outputFilePtr() << force_t_z    << "\t";
            outputFilePtr() << value_f_P_x  << " ";
            outputFilePtr() << value_f_P_y  << " ";
            outputFilePtr() << value_f_P_z  << " ";
            outputFilePtr() << value_f_V_x  << " ";
            outputFilePtr() << value_f_V_y  << " ";
            outputFilePtr() << value_f_V_z  << endl;
        }
    }
    return 0;
}
