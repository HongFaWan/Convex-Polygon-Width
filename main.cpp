#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
//#include "sophus/se3.hpp"


/// 本程序演示sophus的基本用法
//using namespace std;
//using namespace Eigen;
//int main(int argc, char **argv) {
//
//    // 沿Z轴转90度的旋转矩阵
//    Matrix3d R = AngleAxisd(M_PI / 4, Vector3d(0, 0, 1)).toRotationMatrix();
//
//    cout<<R<<endl;
//    // 或者四元数
//    Quaterniond q(R);
//    Sophus::SO3d SO3_R(R);              // Sophus::SO3d可以直接从旋转矩阵构造
//    Sophus::SO3d SO3_q(q);              // 也可以通过四元数构造
//    // 二者是等价的
//    cout << "SO(3) from matrix:\n" << SO3_R.matrix() << endl;
//    cout << "SO(3) from quaternion:\n" << SO3_q.matrix() << endl;
//    cout << "they are equal" << endl;
//
//    // 使用对数映射获得它的李代数
//    Vector3d so3 = SO3_R.log();
//    cout << "so3 = " << so3.transpose() << endl;
//    // hat 为向量到反对称矩阵
//    cout << "so3 hat=\n" << Sophus::SO3d::hat(so3) << endl;
//    // 相对的，vee为反对称到向量
//    cout << "so3 hat vee= " << Sophus::SO3d::vee(Sophus::SO3d::hat(so3)).transpose() << endl;
//
//    // 增量扰动模型的更新
//    Vector3d update_so3(1e0, 0, 0); //假设更新量为这么多
//    Sophus::SO3d SO3_updated = Sophus::SO3d::exp(update_so3) * SO3_R;
//    cout << "SO3 updated = \n" << SO3_updated.matrix() << endl;
//
//    cout << "*******************************" << endl;
//    // 对SE(3)操作大同小异
//    Vector3d t(1, 0, 0);           // 沿X轴平移1
//    Sophus::SE3d SE3_Rt(R, t);           // 从R,t构造SE(3)
//    Sophus::SE3d SE3_qt(q, t);            // 从q,t构造SE(3)
//    cout << "SE3 from R,t= \n" << SE3_Rt.matrix() << endl;
//    cout << "SE3 from q,t= \n" << SE3_qt.matrix() << endl;
//    // 李代数se(3) 是一个六维向量，方便起见先typedef一下
//    typedef Eigen::Matrix<double, 6, 1> Vector6d;
//    Vector6d se3 = SE3_Rt.log();
//    cout << "se3 = " << se3.transpose() << endl;
//    // 观察输出，会发现在Sophus中，se(3)的平移在前，旋转在后.
//    // 同样的，有hat和vee两个算符
//    cout << "se3 hat = \n" << Sophus::SE3d::hat(se3) << endl;
//    cout << "se3 hat vee = " << Sophus::SE3d::vee(Sophus::SE3d::hat(se3)).transpose() << endl;
//
//    // 最后，演示一下更新
//    Vector6d update_se3; //更新量
//    update_se3.setZero();
//    update_se3(0, 0) = 1e-4d;
//
//    Sophus::SE3d SE3_updated = Sophus::SE3d::exp(update_se3) * SE3_Rt;
//    cout << "SE3 updated = " << endl << SE3_updated.matrix() << endl;
//
//    return 0;
//}

Eigen::Vector2d convert(Eigen::Vector2d L,Eigen::Vector3d axis,double angle){
    double r=10;
    Eigen::Matrix3d R=Eigen::AngleAxisd(angle,axis).toRotationMatrix();
    Eigen::Vector3d X;
    Eigen::Vector2d L_;
    X[0]=r*cos(L[1])*cos(L[0]);
    X[1]=r*cos(L[1])*sin(L[0]);
    X[2]=r*sin(L[1]);
    X=R*X;
    L_[0]=atan2(X[1],X[0]);
    L_[1]=atan2(X[2],sqrt(X[1]*X[1]+X[0]*X[0]));

    return L_;
}

double azimuth()

int main(int argc, char **argv){
//    Eigen::Matrix3d R=Eigen::AngleAxisd(M_PI / 4,Eigen::Vector3d(0,0,1)).toRotationMatrix();
    Eigen::Vector3d axis(0,0,1);
    double angle=M_PI/1.0;
    Eigen::Vector2d L(0,M_PI / 4.0);
    Eigen::Vector2d L_;
    L_=convert(L,axis,angle);
    std::cout<<L_;

}

