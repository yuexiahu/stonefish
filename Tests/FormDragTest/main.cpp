/*
    This file is a part of Stonefish.

    Stonefish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Stonefish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

//
//  main.cpp
//  FormDragTest
//
//  无窗口、断言式的形式阻力回归测试。
//  验证完全浸没与跨水面两条阻尼路径的逐面形式阻力，
//  在无水流直线运动时退化为标准阻力方程 0.5*rho*Cd*A_proj*|v|*v。
//  失败时返回非零退出码。
//

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include "core/ConsoleSimulationApp.h"
#include "core/SimulationManager.h"
#include "entities/SolidEntity.h"
#include "entities/forcefields/Ocean.h"
#include "graphics/OpenGLDataStructs.h"

// 无窗口的仿真环境：Entity 构造函数需要通过应用访问名称管理器。
// SimulationManager 只有 BuildScenario 是纯虚函数，场景为空即可。
class HeadlessSimulationManager : public sf::SimulationManager
{
public:
    HeadlessSimulationManager() : SimulationManager(60, sf::Solver::SI, sf::CollisionFilter::EXCLUSIVE) {}
    void BuildScenario() override {}
};

namespace
{

const double RHO = 1000.0;   // 水的密度 [kg/m^3]
const double G_TOL = 1e-2;   // 通用相对容差

int failures = 0;

void check(bool ok, const std::string& what)
{
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", what.c_str());
    if(!ok) ++failures;
}

void checkNear(double got, double expected, double relTol, const std::string& what)
{
    double rel = std::fabs(got - expected) / (std::fabs(expected) > 1e-12 ? std::fabs(expected) : 1.0);
    bool ok = rel <= relTol;
    std::printf("[%s] %s: got %.6f, expected %.6f (rel err %.3e, tol %.1e)\n",
                ok ? "PASS" : "FAIL", what.c_str(), got, expected, rel, relTol);
    if(!ok) ++failures;
}

// 单位立方体（边长 1，中心在原点），每个面两个三角形。
// 绕序保证 cross(p2-p1, p3-p1) 指向外侧。
sf::PlainMesh* MakeUnitCube()
{
    sf::PlainMesh* mesh = new sf::PlainMesh();
    const float h = 0.5f;
    const glm::vec3 pos[8] = {
        glm::vec3(-h,-h,-h), glm::vec3( h,-h,-h), glm::vec3( h, h,-h), glm::vec3(-h, h,-h),
        glm::vec3(-h,-h, h), glm::vec3( h,-h, h), glm::vec3( h, h, h), glm::vec3(-h, h, h)
    };
    for(int i=0; i<8; ++i)
    {
        sf::Vertex v;
        v.pos = pos[i];
        mesh->vertices.push_back(v);
    }

    const unsigned int tri[12][3] = {
        {1,6,5},{1,2,6},   // +X
        {0,7,3},{0,4,7},   // -X
        {3,6,2},{3,7,6},   // +Y
        {0,5,4},{0,1,5},   // -Y
        {4,5,6},{4,6,7},   // +Z
        {0,3,2},{0,2,1}    // -Z
    };
    for(int i=0; i<12; ++i)
    {
        sf::Face f;
        f.vertexID[0] = tri[i][0];
        f.vertexID[1] = tri[i][1];
        f.vertexID[2] = tri[i][2];
        mesh->faces.push_back(f);
    }
    return mesh;
}

// 半径 R 的 UV 球，中心在原点。绕序保证法线指向外侧。
sf::PlainMesh* MakeSphere(float R, int stacks, int slices)
{
    sf::PlainMesh* mesh = new sf::PlainMesh();
    for(int i=0; i<=stacks; ++i)
    {
        float theta = float(M_PI) * float(i) / float(stacks);
        for(int j=0; j<=slices; ++j)
        {
            float phi = 2.f * float(M_PI) * float(j) / float(slices);
            sf::Vertex v;
            v.pos = glm::vec3(R * std::sin(theta) * std::cos(phi),
                              R * std::sin(theta) * std::sin(phi),
                              R * std::cos(theta));
            mesh->vertices.push_back(v);
        }
    }
    auto id = [slices](int i, int j) { return (unsigned int)(i * (slices + 1) + j); };
    for(int i=0; i<stacks; ++i)
    {
        for(int j=0; j<slices; ++j)
        {
            // (i,j) -> (i+1,j) -> (i,j+1) 方向的外法线由参数化导数叉积给出
            sf::Face f1;
            f1.vertexID[0] = id(i, j);
            f1.vertexID[1] = id(i+1, j);
            f1.vertexID[2] = id(i, j+1);
            mesh->faces.push_back(f1);

            sf::Face f2;
            f2.vertexID[0] = id(i+1, j);
            f2.vertexID[1] = id(i+1, j+1);
            f2.vertexID[2] = id(i, j+1);
            mesh->faces.push_back(f2);
        }
    }
    return mesh;
}

// 刚体位于水下，无旋转，沿 +X 平动
sf::Transform MakeBodyTransform()
{
    return sf::Transform(sf::Quaternion::getIdentity(), sf::Vector3(0, 0, 10));
}

// 计算修正后的形式阻力及原始形式阻力矩
void ComputeCorrectedSubmerged(sf::Ocean* ocn, const sf::PlainMesh* mesh, double U,
                               sf::Vector3& Fdq, sf::Vector3& Tdq)
{
    sf::Vector3 Fdf(0,0,0), Tdf(0,0,0);
    sf::SolidEntity::ComputeHydrodynamicForcesSubmerged(mesh, ocn, MakeBodyTransform(), MakeBodyTransform(),
                                                        sf::Vector3(U, 0, 0), sf::V0(), Fdq, Tdq, Fdf, Tdf);
    sf::Vector3 Fdfc(0,0,0), Tdfc(0,0,0);
    sf::SolidEntity::CorrectHydrodynamicForces(ocn, Fdq, Tdq, Fdfc, Tdfc,
                                               sf::Vector3(1,1,1), sf::Vector3(1,1,1), sf::Transform::getIdentity());
}

void ComputeCorrectedSurface(sf::Ocean* ocn, const sf::PlainMesh* mesh, double U,
                             sf::Vector3& Fdq, sf::Vector3& Tdq)
{
    sf::HydrodynamicsSettings settings;
    settings.dampingForces = true;
    settings.reallisticBuoyancy = false;

    sf::Vector3 Fb(0,0,0), Tb(0,0,0), Fdf(0,0,0), Tdf(0,0,0);
    sf::Scalar Swet(0), Vsub(0);
    sf::Renderable debug;
    sf::SolidEntity::ComputeHydrodynamicForcesSurface(settings, mesh, ocn, MakeBodyTransform(), MakeBodyTransform(),
                                                      sf::Vector3(U, 0, 0), sf::V0(),
                                                      Fb, Tb, Fdq, Tdq, Fdf, Tdf, Swet, Vsub, debug);
    sf::Vector3 Fdfc(0,0,0), Tdfc(0,0,0);
    sf::SolidEntity::CorrectHydrodynamicForces(ocn, Fdq, Tdq, Fdfc, Tdfc,
                                               sf::Vector3(1,1,1), sf::Vector3(1,1,1), sf::Transform::getIdentity());
}

double Magnitude(const sf::Vector3& v)
{
    return v.length();
}

bool OppositeToMotion(const sf::Vector3& F)
{
    double m = Magnitude(F);
    return m > 0 && F.getX() < 0 && std::fabs(F.getY()) < 1e-3 * m && std::fabs(F.getZ()) < 1e-3 * m;
}

double RawForceMagnitudeSubmerged(sf::Ocean* ocn, const sf::PlainMesh* mesh, double U)
{
    sf::Vector3 Fdq(0,0,0), Tdq(0,0,0), Fdf(0,0,0), Tdf(0,0,0);
    sf::SolidEntity::ComputeHydrodynamicForcesSubmerged(mesh, ocn, MakeBodyTransform(), MakeBodyTransform(),
                                                        sf::Vector3(U, 0, 0), sf::V0(), Fdq, Tdq, Fdf, Tdf);
    return Magnitude(Fdq);
}

double RawTorqueMagnitudeSubmerged(sf::Ocean* ocn, const sf::PlainMesh* mesh, double U)
{
    sf::Vector3 Fdq(0,0,0), Tdq(0,0,0), Fdf(0,0,0), Tdf(0,0,0);
    sf::SolidEntity::ComputeHydrodynamicForcesSubmerged(mesh, ocn, MakeBodyTransform(), MakeBodyTransform(),
                                                        sf::Vector3(U, 0, 0), sf::V0(), Fdq, Tdq, Fdf, Tdf);
    return Magnitude(Tdq);
}

} // namespace

int main(int argc, const char* argv[])
{
    (void)argc; (void)argv;

    HeadlessSimulationManager manager;
    sf::ConsoleSimulationApp app("FormDragTest", DATA_DIR_PATH, &manager); // 无窗口，不启动仿真循环

    sf::Fluid water;
    water.name = "Water";
    water.density = RHO;
    water.viscosity = 1e-6;
    water.IOR = 1.33;
    sf::Ocean ocn("FormDragTestOcean", 0.0, water); // 无波浪、无海流

    sf::PlainMesh* cube = MakeUnitCube();
    sf::PlainMesh* sphere = MakeSphere(1.f, 32, 64);

    std::printf("== 单位立方体：标准阻力方程 F = 0.5*rho*Cd*A_proj*U^2 ==\n");
    // A_proj = 1 m^2, Cd = 1, rho = 1000 -> F = 500*U^2
    const double cubeExpected[3] = {125.0, 500.0, 2000.0};
    const double cubeU[3] = {0.5, 1.0, 2.0};
    double cubeF[3];
    for(int i=0; i<3; ++i)
    {
        sf::Vector3 F, T;
        ComputeCorrectedSubmerged(&ocn, cube, cubeU[i], F, T);
        cubeF[i] = Magnitude(F);
        checkNear(cubeF[i], cubeExpected[i], G_TOL,
                  "submerged cube U=" + std::to_string(cubeU[i]) + " m/s force [N]");
        check(OppositeToMotion(F),
              "submerged cube U=" + std::to_string(cubeU[i]) + " m/s force opposes motion");
    }

    std::printf("== 速度缩放 ==\n");
    checkNear(cubeF[2] / cubeF[1], 4.0, G_TOL, "F(2U)/F(U)");
    checkNear(cubeF[0] / cubeF[1], 0.25, G_TOL, "F(0.5U)/F(U)");

    std::printf("== 对称立方体绕质心的原始阻力矩应可忽略 ==\n");
    {
        double Fmag = RawForceMagnitudeSubmerged(&ocn, cube, 1.0);
        double Tmag = RawTorqueMagnitudeSubmerged(&ocn, cube, 1.0);
        check(Tmag < 1e-3 * Fmag * 0.5 * 10.0,
              "submerged cube torque about CG negligible (|T|=" + std::to_string(Tmag) + ")");
    }

    std::printf("== 半径 1 m 球体，U=2 m/s ==\n");
    {
        sf::Vector3 F, T;
        ComputeCorrectedSubmerged(&ocn, sphere, 2.0, F, T);
        double expected = 0.5 * RHO * 1.0 * M_PI * 1.0 * 4.0; // 0.5*rho*Cd*pi*R^2*U^2
        checkNear(Magnitude(F), expected, 0.01, "submerged sphere force [N]");
        check(OppositeToMotion(F), "submerged sphere force opposes motion");
    }

    std::printf("== 全浸没立方体：surface 与 submerged 两条路径一致 ==\n");
    for(int i=0; i<3; ++i)
    {
        sf::Vector3 Fs, Ts, Fc, Tc;
        ComputeCorrectedSurface(&ocn, cube, cubeU[i], Fs, Ts);
        ComputeCorrectedSubmerged(&ocn, cube, cubeU[i], Fc, Tc);
        double diff = Magnitude(Fs - Fc);
        double ref = Magnitude(Fc) > 1e-12 ? Magnitude(Fc) : 1.0;
        check(diff / ref < 1e-3,
              "U=" + std::to_string(cubeU[i]) + " m/s surface/submerged match (|dF|=" + std::to_string(diff) + ")");
        check(OppositeToMotion(Fs), "U=" + std::to_string(cubeU[i]) + " m/s surface force opposes motion");
    }

    delete cube;
    delete sphere;

    if(failures > 0)
    {
        std::printf("\n%d assertion(s) FAILED\n", failures);
        return 1;
    }
    std::printf("\nAll assertions passed\n");
    return 0;
}
