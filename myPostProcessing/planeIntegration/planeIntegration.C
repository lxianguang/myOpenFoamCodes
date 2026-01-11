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
    planeIntegration

Group
    postProcessing

Description
    截取涡量场yz平面并进行面积分计算

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "fvCFD.H"
#include "vector.H"
#include "argList.H"
#include "OFstream.H"
#include "polyPatchID.H"
#include "sampledPlane.H"
#include "surfaceFields.H"
#include "volPointInterpolation.H" 

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // 准备参数列表
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"   

    Info << "================ Loading integration parameters ===============" << endl;
    Info << "Loading parameters from " << runTime.system() + "/integrationControlDict" << endl;

    IOdictionary integrationControlDict
    (
        IOobject
        (
            "integrationControlDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    const scalar myruntime   = integrationControlDict.get<scalar>("runTime");   // 选取的积分时刻
    const scalar startX      = integrationControlDict.get<scalar>("startX");    // 积分起始位置
    const scalar interval    = integrationControlDict.get<scalar>("interval");  // 数据采样间隔
    const scalar endX        = integrationControlDict.get<scalar>("endX");      // 积分结束位置
    const scalar raidus      = integrationControlDict.get<scalar>("radius");    // 积分区域半径
    const scalar tolerance   = integrationControlDict.get<scalar>("tolerance"); 
    const fileName outputDir = "./postProcessing/integration/"; // 定义积分时间与输出文件路径

    instantList timeDirs(1);
    timeDirs[0] = instant(myruntime);

    // 创建输出文件及写入文件头
    if (!isDir(outputDir)) mkDir(outputDir);
    autoPtr<OFstream> outputFilePtr;
    outputFilePtr.reset(new OFstream(outputDir/"planeIntegration.dat"));
    outputFilePtr() << "Variables = coordinate, area, indicator1, indicator2, indicator3, indicator4, indicator4, omega_x, omega_y, omega_z, drag, lift" << "\n" << endl;

    Info << "==================== Loading velocity field ===================" << endl;
    forAll(timeDirs, timeI)
    {
        mesh.readUpdate();                       // 更新网格信息
        runTime.setTime(timeDirs[timeI], timeI); // 读取文件时刻

        volVectorField velocity                  // 读取速度场
        (
            IOobject
            (
                "U",
                runTime.timeName(), 
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );                

        Info << "================= Calculating vorticity field =================" << endl;
        volVectorField omega = fvc::curl(velocity);             // 计算涡量场

        Info << "====================== Cutting yz planes ======================" << endl;
        label nPlanes = label((endX - startX) / interval) + 1;  // 根据interval参数创建多个平面
        for (label planeI = 0; planeI < nPlanes; planeI++)
        {
            scalar xPos = startX + planeI * interval;           // 计算当前平面的x坐标
            point  planePoint(xPos, 0, 0);
            vector normal(1, 0, 0);                             // 法向量指向x方向
            
            Info << "cutting YZ plane [" << planeI + 1 << "/" << nPlanes << "]: x = " << xPos << " -----------------------------------" << endl;

            sampledPlane yzPlane // 创建采样平面
            (
                "yzPlane",
                mesh,
                plane(planePoint, normal)
            );
            yzPlane.update();

            // 创建插值对象
            autoPtr<interpolation<vector>> interpOmegaPtr(
                interpolation<vector>::New("cellPoint", omega)
            );

            const vectorField omegaSampled = yzPlane.sample(interpOmegaPtr())();  // 使用sample方法获取平面上的涡量场
            const scalarField areasSampled = yzPlane.magSf();                     // 获取每个面片的面积
            const vectorField planeCenters = yzPlane.Cf();                        // 获取面片中心坐标

            // 进行面积分计算
            scalar dragforces = 0.0;
            scalar liftforces = 0.0;
            scalar totalAreas = 0.0;
            scalar summation1 = 0.0;
            scalar summation2 = 0.0;
            scalar summation3 = 0.0;
            scalar summation4 = 0.0;
            scalar summation5 = 0.0;

            vector averOmegas = sum(omegaSampled * areasSampled) / sum(areasSampled);

            forAll(omegaSampled, faceI)
            {
                dragforces += -500 * (planeCenters[faceI].y() * omegaSampled[faceI].z() - planeCenters[faceI].z() * omegaSampled[faceI].y()) * areasSampled[faceI];
                liftforces += -1000 * planeCenters[faceI].z() * omegaSampled[faceI].x() * areasSampled[faceI];
                if (mag(planeCenters[faceI] - point(xPos, 0, 0)) <= raidus)
                {
                    summation1 += mag(omegaSampled[faceI].x()) * areasSampled[faceI];
                    summation2 += mag(omegaSampled[faceI].y()) * areasSampled[faceI];
                    summation3 += mag(omegaSampled[faceI].z()) * areasSampled[faceI];
                    summation4 += mag(omegaSampled[faceI].z() - averOmegas.z()) * areasSampled[faceI];
                    if (mag(omegaSampled[faceI].y()) > tolerance && mag(omegaSampled[faceI].z()) > tolerance)
                    {
                        summation5 += mag(omegaSampled[faceI].x()/(omegaSampled[faceI].y())) * areasSampled[faceI];
                        totalAreas += areasSampled[faceI];
                    }
                }
            }

            scalar indicator1 = summation5 / totalAreas;
            scalar indicator2 = summation4 / summation2;
            scalar indicator3 = summation3 / summation2;
            scalar indicator4 = summation1 / summation2;
            
            Info << "absulote summations : " << summation1 << "    " << summation2 << "    " << summation3 << "    " << summation4 << "    " << summation5 << endl;
            Info << "area and outputs    : " << totalAreas << "    " << indicator1 << "    " << indicator2 << "    " << indicator3 << "    " << indicator4 << endl;
            
            // 输出结果
            outputFilePtr() << xPos << "    "
                            << totalAreas << "    "
                            << indicator1 << "    "
                            << indicator2 << "    "
                            << indicator3 << "    "
                            << indicator4 << "    "
                            << summation1 << "    "
                            << summation2 << "    "
                            << summation3 << "    "
                            << dragforces << "    "
                            << liftforces << endl;
        }
    }

    return 0;
}
