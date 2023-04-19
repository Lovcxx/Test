close all
clear
clc
disp('***** 基于EKF的位置速度观测组合导航程序 *****');
disp('Step1:加载数据;');

load IMU_data200.mat       %惯导原始数据
load Reference_data.mat    %GPS测量数据

disp('Step2:初始化参数;');
%% 一些导航参数常数项
WIE = 7.292115e-5;           % 地球自转角速度
r0  = 6378137.0;             % 地球半径
EE  =  0.0818191908426;      % 偏心率
d2r   =  pi/180;             % degree to radian角度转弧度
r2d   =  180/pi;             % radian to degree弧度转角度
dh2rs =  d2r/3600;           % deg/h to rad/s角速度转弧速度
%% 导航坐标系下初始化姿态，速度，位置
yaw = (0)*pi/180;%航向角
pitch = 0*pi/180;%俯仰角
roll = 0*pi/180;%滚动角
cbn=eul2dcm(roll,pitch,yaw);
cnb=cbn';
q=dcm2quat(cbn)';

Vn=0;%北向速度
Ve=0;%东向速度
Vd=0;%地向速度
V_last=[Vn Ve Vd]';

Lati = 31.4913627505302*pi/180;%纬度
Longi= 120.849577188492*pi/180;%经度
Alti = 6.6356;%高度

sampt0=1/200;%惯导系统更新时间

Rn = r0*(1-EE^2)/(1-EE^2*(sin(Lati))^2)^1.5;         %子午圈曲率半径
Re = r0/(1-EE^2*(sin(Lati))^2)^0.5;                  %卯酉圈曲率半径
g_u = -9.7803267711905*(1+0.00193185138639*sin(Lati)^2)...
    /((1-0.00669437999013*sin(Lati)^2)^0.5 *(1.0 + Alti/r0)^2);
g = [0 0 -g_u]';%重力
g0=9.80665;
%% 卡尔曼滤波P、Q、R设置
% P的设置（估计的协方差矩阵）
std_roll    =   (5)*d2r;
std_pitch   =   (5)*d2r;
std_yaw     =   (60)*d2r;

std_vel     =   0.1;
std_pos     =   5;

std_gyro    =   3*0.5*dh2rs;             % 陀螺随机漂移0.5度/小时
std_acc     =   3*0.15e-3*g0;            % 加表零偏0.15mg

Pfilter     =   diag([std_roll^2 std_pitch^2 std_yaw^2 std_vel^2 std_vel^2 std_vel^2 (std_pos/3600/30/57.3)^2 (std_pos/3600/30/57.3)^2 std_pos^2  std_gyro^2 std_gyro^2 std_gyro^2 std_acc^2 std_acc^2 std_acc^2]);

% Q的设置（过程噪声的协方差矩阵）
std_Wg      =   0.15*(2.909*1e-4);       % 陀螺漂移噪声,度/根号小时转化成rad/根号秒
std_Wa      =   0.21/60/3;               % 加表漂移噪声
Qkf         =   diag([std_Wg^2 std_Wg^2 std_Wg^2 std_Wa^2 std_Wa^2 std_Wa^2]);

G			=   zeros(15, 6);
F = zeros(15);
F_i=zeros(9,9);
F_s=zeros(9,6);
H           =   zeros(6,15);
H(1:3,4:6)  =   eye(3);
H(4:6,7:9)  =   eye(3);
% R的设置（量测噪声协方差矩阵）
R           =   diag([std_vel^2 std_vel^2 std_vel^2 (std_pos/3600/30/57.3)^2 (std_pos/3600/30/57.3)^2 (std_pos)^2]);
%初始化一些参数
cnt_GPS=1;                          % GPS计数器
k_StartGPS=15;                      % GPS开始计算位置的帧数
X_kalman=zeros(size(Pfilter,1),1);  % 状态量
P_kalman=Pfilter;                   % 构造状态协方差矩阵
output_data=cell(1,length(IMU));    % 保存EKF矫正后的导航信息

% 开始计算
for cnt_IMU=1:length(IMU) %循环运行程序
    if mod(cnt_IMU,5000)==0 %间隔5000次输出运行进度
       disp(['Time: ', num2str(IMU{cnt_IMU}.time)]);
    end
    
    % 获取IMU数据
    WM=IMU{cnt_IMU}.gyro'-X_kalman(10:12);%去除零偏
    AM=IMU{cnt_IMU}.acce'-X_kalman(13:15);%计算真实测量值
    
    % 动力学状态方程的状态转移矩阵F，连续型
    F_i(1,2)=WM(3);F_i(1,3)=-WM(2);
    F_i(2,1)=-WM(3);F_i(2,3)=WM(1);
    F_i(3,1)=WM(2);F_i(3,2)=-WM(1);
    F_i(1:3,4:6)=cbn;
    F_s(1,1)=1;F_s(2,2)=1;F_s(3,3)=1;
    F_s(4,1)=g(1);F_s(4,2)=g(2);F_s(4,3)=g(3);
    F_s(5,1)=V_last(2)/(Re+Alti);F_s(5,2)=-V_last(1)/(Rn+Alti);F_s(5,3)=-V_last(2)*tan(Lati)/(Re+Alti);
    F_s(6,1)=V_last(1)/((Rn+Alti)*cos(Lati));F_s(6,2)=V_last(2)/((Re+Alti)*cos(Lati));F_s(6,3)=-V_last(1)*tan(Lati)/((Rn+Alti)*cos(Lati));
    F_s(7,4)=1;F_s(8,5)=1;F_s(9,6)=1;
    F=[eye(9)+F_i*sampt0, F_s*sampt0, zeros(9,1); zeros(1,9),1];
    % 动力学状态转移矩阵的协方差矩阵Q
    Q_cont      =   Qkf*sampt0;
    Q_kalman =   [[Q_cont, zeros(6,3)];[zeros(3,6), zeros(3,3)]];
    
    % 更新状态量
    X_kalman=F*X_kalman+[zeros(9,1);g];
    % 计算姿态让update函数使用
    [roll, pitch, yaw]=dcm2angle(cbn,'zyx');
    % 低精度姿态解算
    % 角速度积分得到姿态
    [cbn, q]=gyro_update(cbn, WM, sampt0, q);
    % 加速度向量的低通滤波得到姿态
    [cbn, q]=accel_update(cbn, AM, sampt0, q, roll, pitch);
    % 解算卡尔曼滤波的状态估计及方差估计    
    P_kalman=F*P_kalman*F'+Q_kalman; %状态预测协方差矩阵
    K=P_kalman*H'*inv(H*P_kalman*H'+R);%计算卡尔曼增益矩阵
    % 更新状态变量
    if cnt_GPS<=length(Reference_data) && IMU{cnt_IMU}.time==Reference_data{cnt_GPS}.time
        % 如果当前时刻有GPS观测值
        [vn,ve,vr,LL]=BLH2xyzV(Reference_data{cnt_GPS}.latitude*pi/180,Reference_data{cnt_GPS}.longitude*pi/180,Reference_data{cnt_GPS}.altitude,0,0,0,0,0,0,0,0,jd);%给出参考位置
        z_gps=[vn-IMU{cnt_IMU}.vn; ve-IMU{cnt_IMU}.ve; vr-IMU{cnt_IMU}.vd;LL(1)*r2d*3600;LL(2)*r2d*3600;LL(3)];
        X_kalman=X_kalman+K*(z_gps(1:6)-H(1:6,:)*X_kalman);%状态变量更新
        P_kalman=(eye(size(K,1))-K*H)*P_kalman;%协方差矩阵更新
        cnt_GPS=cnt_GPS+1;%GPS计数器自增
    else
        % 如果当前时刻没有GPS观测值
        X_kalman=X_kalman;%状态变量更新
        P_kalman=P_kalman;%协方差矩阵更新
    end
    
    %保存所有导航信息
    output_data{cnt_IMU}.time=IMU{cnt_IMU}.time;
    output_data{cnt_IMU}.roll=roll*r2d;
    output_data{cnt_IMU}.pitch=pitch*r2d;
    output_data{cnt_IMU}.yaw=yaw*r2d;
    output_data{cnt_IMU}.Vn=X_kalman(4);
    output_data{cnt_IMU}.Ve=X_kalman(5);
    output_data{cnt_IMU}.Vd=X_kalman(6);
    output_data{cnt_IMU}.lati=LL(1);
    output_data{cnt_IMU}.longi=LL(2);
    output_data{cnt_IMU}.alti=LL(3);
    V_last=[X_kalman(4),X_kalman(5),X_kalman(6)];
end