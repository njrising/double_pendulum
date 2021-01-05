#include "OpenGL.h"
// Simulation of double pendulum
// based on modeling techniques from paper
//"Physically Based Modeling: Principles and Practice" by Andrew Witkin
// http://www.cs.cmu.edu/~baraff/sigcourse/notesf.pdf
// Author: Nathan Rising
// Last Revision: 7/20/2020

//----------------------------- Begin simulation ------------------------------/
const float length = 0.5;                              // Size of rigid links
const float c = cos(-0.5) * length;                    // Start cosine
const float s = sin(-0.5) * length;                    // Start sine
float globalPosition[4]  = {c,s, 2 * c,2 * s};         // Particles position
float globalVelocity[4] = {0.0,0.0,0.0,0.0};           // Particles velocity
float globalAcceleration[4] = {0.0,0.0,0.0,0.0};       // Particles acceleration
//float globalMass[2][2];                              // Particles Mass (2X2)
float fExt[] = {0.0,-0.5, 0.0,-0.5};                   // External force
float dt = 0.04;                                       // Time step increment
vec2 pos1 = {&globalPosition[0], &globalPosition[1]};
vec2 pos2 = {&globalPosition[2], &globalPosition[3]};
vec2 vel1 = {&globalVelocity[0], &globalVelocity[1]};
vec2 vel2 = {&globalVelocity[2], &globalVelocity[3]};
std::vector<float> C = {0.0,0.0};
std::vector<float> cDot(2);
std::vector<std::vector<float>> J(2, std::vector<float>(4, 0.0));
std::vector<std::vector<float>> jDot(2, std::vector<float>(4, 0.0));

std::vector<float> ForceConstraint(){

    J[0] = {2 * pos1.x(), 2 * pos1.y(), 0.0, 0.0};         // (2X4)
    J[1] = {2 * (pos1.x() - pos2.x()), 2 * (pos1.y() - pos2.y()), 2 * (pos2.x() - pos1.x()), 2 * (pos2.y() - pos1.y())};

    jDot[0] = {2 * vel1.x(), 2 * vel2.y(), 0.0, 0.0};      // (2X4)
    jDot[1] = {2 * (vel1.x() - vel2.x()), 2 * (vel1.y() - vel2.y()), 2 * (vel2.x() - vel1.x()), 2 * (vel2.y() - vel1.y())};

    C[0] = dot(pos1,pos1) - length * length;               // (2X1)
    C[1] = dot(pos1,pos1) + dot(pos2,pos2) - 2 * dot(pos1,pos2) - length * length;

    cDot[0] = 2.0 * dot(pos1, vel1);                       // (4X1)
    cDot[1] = 2.0 * (dot(pos1,vel1) + dot(pos2,vel2) - dot(pos1,vel2) - dot(pos2,vel1));

    float ks = 300.0;                                      // spring constant
    float kd = 1.0;                                        // damping constant

    // A = J * W * J'
    std::vector<std::vector<float>> A(2,std::vector<float>(2));
    for(int i = 0;i<2;++i){
        for(int j = 0;j<2;++j){
            for(int k = 0;k<4;++k){
                A[i][j] += J[i][k] * J[j][k];}}}

    // B = -Jdot * globalVelocity - J * W * fExt - ks * Constraint - kd * ConstraintDot
    std::vector<float> B(2, 0.0);
    for(unsigned k = 0;k < 2;++k){
        B[k] += -C[k] * ks - cDot[k] * kd;
        for(unsigned l = 0;l < 4;++l){
        B[k] += -jDot[k][l] * globalVelocity[l];     // B = jDot * globalVelocity
        B[k] += -J[k][l] * fExt[l];}}                // J * W * fExt

    //solving Ax = B
    std::vector<float> Xo = {0.0,0.0};

    for(int i = 0;i < 40;++i){
        Xo[0] = -1.0/A[0][0] * (A[0][1] * Xo[1] - B[0]);
        Xo[1] = -1.0/A[1][1] * (A[1][0] * Xo[0] - B[1]);
    }

    // J' * Xo
    std::vector<float> FConstraint(6);
    FConstraint[0] = J[0][0] * Xo[0] + J[1][0] * Xo[1];
    FConstraint[1] = J[0][1] * Xo[0] + J[1][1] * Xo[1];
    FConstraint[2] = J[0][2] * Xo[0] + J[1][2] * Xo[1];
    FConstraint[3] = J[0][3] * Xo[0] + J[1][3] * Xo[1];

return FConstraint;
}

int main(){

OpenGL program;

while(!program.close()){
    std::vector<float> Fc = ForceConstraint();
    for(int i = 0;i < 4;++i){
        globalAcceleration[i] = fExt[i] + Fc[i];
        globalVelocity[i] += globalAcceleration[i] * dt;
        globalPosition[i] += globalVelocity[i] * dt;
    }
    /* Limit rendering frames */
    usleep(10000);
    /* Render */
    program.poll();
    program.render({0.0,0.0,globalPosition[0],globalPosition[1],globalPosition[2],globalPosition[3]});
    program.post();
}
glfwTerminate();
return 0;
}
