#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<algorithm>
#include <fstream>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <algorithm>
#include <iomanip>
using namespace std;
struct node
{
    double x,y;
}p[90000],sta[90000];

int j,n,top;
double js;
double Radius=1.0;

double multi(node p1,node p2, node p0)
{
    double x1=p1.x-p0.x;
    double y1=p1.y-p0.y;

    double x2=p2.x-p0.x;
    double y2=p2.y-p0.y;

    return(x1*y2-x2*y1);
}//计算叉积

double dis(node p1,node p2)
{
    return (double((p1.x-p2.x)*(p1.x-p2.x))+double((p1.y-p2.y)*(p1.y-p2.y)));
}//计算任意两点间距离


double Vertical_dis(node pi1,node pi2,node pj)
{
    cout<<sqrt((pi1.x-pi2.x)*(pi1.x-pi2.x)+(pi1.y-pi2.y)*(pi1.y-pi2.y))<<endl;
    return fabs(multi(pi1, pi2, pj))/sqrt((pi1.x-pi2.x)*(pi1.x-pi2.x)+(pi1.y-pi2.y)*(pi1.y-pi2.y));
}//计算点j到直线ii和i2的垂直距离

bool cmp(node p1,node p2)
{
    double tt=multi(p1,p2,p[1]);
    if(tt<0)return(0);
    if(tt==0 && dis(p1,p[1])>dis(p2,p[1]))return(0);
    return(1);
}

void graham()
{
    sort(p+2,p+n+1,cmp);
    for(int i=1;i<=2;i++)sta[i]=p[i];
    top=2;
    n++;p[n]=p[1];

    for(int i=3;i<=n;i++)
    {
        while(top>1 && multi(sta[top],p[i],sta[top-1])<=0)top--;
        top++;
        sta[top]=p[i];

    }
}//寻找凸多边形的顶点

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
}//

Eigen::Vector3d compute(node focus,double azimuth){
    Eigen::Matrix3d R1=Eigen::AngleAxisd(-focus.x,Eigen::Vector3d(0,0,1)).toRotationMatrix();
    Eigen::Matrix3d R2=Eigen::AngleAxisd(focus.y,Eigen::Vector3d(0,1,0)).toRotationMatrix();
    Eigen::Matrix3d R3=Eigen::AngleAxisd(M_PI/2.0,Eigen::Vector3d(0,0,1)).toRotationMatrix();
    Eigen::Matrix3d R4=Eigen::AngleAxisd(-azimuth,Eigen::Vector3d(0,1,0)).toRotationMatrix();
    Eigen::Matrix3d R=R4*R3*R2*R1;

    cout<<R<<endl;
    double y=acos(R(2,2));
//    double z2_1=acos(R(2,0)/(-sin(y)));
//    double z2_2=asin(R(2,1)/(sin(y)));
    double z1=atan2(R(2,1),-R(2,0));
    double z2=atan2(R(1,2),R(0,2));

    Eigen::Matrix3d Rz1=Eigen::AngleAxisd(z1,Eigen::Vector3d(0,0,1)).toRotationMatrix();
    Eigen::Matrix3d Rz2=Eigen::AngleAxisd(z2,Eigen::Vector3d(0,0,1)).toRotationMatrix();
    Eigen::Matrix3d Ry=Eigen::AngleAxisd(y,Eigen::Vector3d(0,1,0)).toRotationMatrix();

    Eigen::Vector3d ans(z1,y,z2);
    cout<<Rz2*Ry*Rz1<<endl;
//    cout<<"z2_1="<<z2_1<<"    "<<"z2_2="<<z2_2<<"    "<<"z2_3="<<z2_3<<endl;

    cout<<"azimuth = "<<azimuth<<endl;
    return ans;

}

node project(node p){
    node p_=p;
    return p_;
}//投影方法

double azimuth(node p1,node p2){
    double dlon=fabs(p2.x-p1.x)*M_PI/180.0;
    double dlat=(p2.y-p1.y)*M_PI/180.0;
    double a=sin(dlat/2.0)*sin(dlat/2.0)+cos(p1.y*M_PI/180.0)*cos(p2.y*M_PI/180.0)*sin(dlon/2.0)*sin(dlon/2.0);
    double c=fabs(2*Radius*atan2(sqrt(a),sqrt(1-a)));

    double l_lat=Radius*dlat;
    if(p2.x>=p1.x){
        return asin(sin(l_lat)/sin(c));
    }
    else{
        return -asin(sin(l_lat)/sin(c));
    }
}
//double azimuth(node p1,node p2){
//    double dlon=fabs(p2.x-p1.x);
//    double dlat=(p2.y-p1.y);
//    double a=sin(dlat/2.0)*sin(dlat/2.0)+cos(p1.y)*cos(p2.y)*sin(dlon/2.0)*sin(dlon/2.0);
//    double c=fabs(2*Radius*atan2(sqrt(a),sqrt(1-a)));
//
//    double l_lat=Radius*dlat;
//    if(p2.x>=p1.x){
//        return asin(sin(l_lat)/sin(c));
//    }
//    else{
//        return -asin(sin(l_lat)/sin(c));
//    }
//}

//double azimuth(node p1,node p2){
//    double lat1=p1.y*M_PI/180.0;
//    double lat2=p2.y*M_PI/180.0;
//    double dlon=(p2.x-p1.x)*M_PI/180.0;
//    double dlat=(p2.y-p1.y)*M_PI/180.0;
//
////    return atan2(sin(dlon)*cos(lat2),cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dlon));
//    return atan2(cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dlon),sin(dlon)*cos(lat2));
//
//}

Eigen::Vector3d rotation(string infile_path,string outfile_path){
    ifstream file(infile_path);
    string s;
    vector<Eigen::Vector2d,Eigen::aligned_allocator<Eigen::Vector2d>> X;
    Eigen::Vector2d x;
    vector<Eigen::Vector2d,Eigen::aligned_allocator<Eigen::Vector2d>> Index;
    vector<Eigen::Vector2d,Eigen::aligned_allocator<Eigen::Vector2d>>::iterator it;
    vector<int> len;
    vector<int>::iterator itd;
    node focus;
    focus.x=0.0;
    focus.y=0.0;
    if (!file.is_open()) {
        cout << "Error   opening file";
        exit(1);
    }
    while (!file.eof()) {
        file>>std::fixed>>setprecision(30) >> x[1];
        file>>std::fixed>>setprecision(30) >> x[0];

//        x[0]*=M_PI/180.0;
//        x[1]*=M_PI/180.0;
        X.push_back(x);
    }
    int k=1;
    for(it=X.begin();it!=X.end();it++){
        p[k].x=(*it)[0];
        p[k].y=(*it)[1];
        k++;
    }
    n=X.size();
    for (int i = 1; i <= n; i++) {
        if (p[i].y < p[1].y || (p[i].y == p[1].y && p[i].x < p[1].x)) {
            node t = p[1];
            p[1] = p[i];
            p[i] = t;
        }//找最下方的点

        

    }
    for(int i=1;i<=n;i++){
        focus.x+=p[i].x;
        focus.y+=p[i].y;
    }
    focus.x/=n;
    focus.y/=n;

    cout<<"focus:("<<focus.x<<","<<focus.y<<")"<<endl;
    focus.x*=M_PI/180.0;
    focus.y*=M_PI/180.0;


    top = 0;//记录凸包上点的个数
    graham();

    ofstream outfile(outfile_path);
    for(int i=1;i<=top;i++){
        outfile<<sta[i].x<<"  "<<sta[i].y<<endl;
    }
    outfile.close();


    sta[top + 1] = sta[1];
    j = 2;
    for (int i = 1; i <= top; i++) {
        while (fabs(multi(sta[i], sta[i + 1], sta[j])) < fabs(multi(sta[i], sta[i + 1], sta[j + 1]))) {
            j++;
            if (j >= top)j = 1;
        }
//        js = Vertical_dis(sta[i], sta[i+1],sta[j]);
        js=dis(sta[i],sta[j]);

        len.push_back(js);
        Eigen::Vector2d index;
        index[0]=i;index[1]=j;
        Index.push_back(index);
//        if (js > ans){
//            ans = js;
//            m=i;
//        }
    }//旋转卡壳

    auto element=min_element(len.begin(),len.end());
    int index=distance(len.begin(),element);
    int i1=Index[index][0];
    int i2=Index[index][0]+1;
    int j=Index[index][1];

    cout << "min_width = "<<len[index]<<"  index_i="<<Index[index][0]<<"  index_j="<<Index[index][1]<< endl;
    cout << "i="<<i1<<","<< "i+1="<<i2<< endl;

    cout<<i1<<" point:("<<sta[i1].x<<","<<sta[i1].y<<")"<<endl;
    cout<<i2<<" point:("<<sta[i2].x<<","<<sta[i2].y<<")"<<endl;
    cout<<j<<" point:("<<sta[j].x<<","<<sta[j].y<<")"<<endl;




    cout<<"rotation_angle="<<atan2((sta[i1].y-sta[i2].y),(sta[i1].x-sta[i2].x))*180.0/M_PI<<endl;
    Eigen::Vector3d ans=compute(focus,azimuth(sta[i1],sta[i2]));
    cout<<"rotate_order : z-y-z     after=[2]*[1]*[0]*coordinate"<<endl;
    cout<<setprecision(9) << std::fixed<<ans<<endl;

}

int main() {
//    ifstream file("/home/whf/Macro/sphere_convert/coor2.txt");
//    string s;
//    vector<Eigen::Vector2d,Eigen::aligned_allocator<Eigen::Vector2d>> X;
//    Eigen::Vector2d x;
//    vector<Eigen::Vector2d,Eigen::aligned_allocator<Eigen::Vector2d>> Index;
//    vector<Eigen::Vector2d,Eigen::aligned_allocator<Eigen::Vector2d>>::iterator it;
//    vector<int> len;
//    vector<int>::iterator itd;
//    node focus;
//    focus.x=0;
//    focus.y=0;
//    if (!file.is_open()) {
//        cout << "Error   opening file";
//        exit(1);
//    }
//    while (!file.eof()) {
//        file>>x[0];
//        file>>x[1];
//        X.push_back(x);
//    }
//    int k=1;
//    for(it=X.begin();it!=X.end();it++){
//        p[k].x=(*it)[0];
//        p[k].y=(*it)[1];
//        k++;
//    }
//    n=X.size();
//    for (int i = 1; i <= n; i++) {
//        if (p[i].y < p[1].y || (p[i].y == p[1].y && p[i].x < p[1].x)) {
//            node t = p[1];
//            p[1] = p[i];
//            p[i] = t;
//        }//找最下方的点
//
//        focus.x+=p[i].x;
//        focus.y+=p[i].y;
//
//    }
//    focus.x/=n;
//    focus.y/=n;
//
//    top = 0;//记录凸包上点的个数
//    graham();
//
//    ofstream outfile("/home/whf/Macro/sphere_convert/out2.txt");
//    for(int i=1;i<=top;i++){
//        outfile<<sta[i].x<<"  "<<sta[i].y<<endl;
//    }
//    outfile.close();
//
//
//    sta[top + 1] = sta[1];
//    j = 2;
//    for (int i = 1; i <= top; i++) {
//        while (fabs(multi(sta[i], sta[i + 1], sta[j])) < fabs(multi(sta[i], sta[i + 1], sta[j + 1]))) {
//            j++;
//            if (j >= top)j = 1;
//        }
//        js = dis(sta[i], sta[j]);
//
//        len.push_back(js);
//        Eigen::Vector2d index;
//        index[0]=i;index[1]=j;
//        Index.push_back(index);
//    }//旋转卡壳
//
//    auto element=min_element(len.begin(),len.end());
//    int index=distance(len.begin(),element);
//    int i1=Index[index][0];
//    int i2=Index[index][0]+1;
//
//    cout << len[index]<<"  "<<Index[index][0]<<"  "<<Index[index][1]<< endl;
//    cout << "i="<<i1<<","<< "i+1="<<i2<< endl;
//
//    cout<<Index[index][0]<<" point:("<<sta[int(Index[index][0])].x<<","<<sta[int(Index[index][0])].y<<")"<<endl;
//    cout<<Index[index][1]<<" point:("<<sta[int(Index[index][1])].x<<","<<sta[int(Index[index][1])].y<<")"<<endl;
//
//
//
//
//    cout<<"ans="<<atan2((sta[i1].y-sta[i2].y),(sta[i1].x-sta[i2].x))*180.0/M_PI<<endl;
//    compute(focus,azimuth(sta[i1],sta[i2]));


    rotation("/home/whf/Macro/sphere_convert/Japan_BL.dat","/home/whf/Macro/sphere_convert/Japan_BL_out.txt");



    return 0;
}