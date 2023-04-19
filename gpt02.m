close all
clear
clc
disp('***** ����EKF��λ���ٶȹ۲���ϵ������� *****');
disp('Step1:��������;');

load IMU_data200.mat       %�ߵ�ԭʼ����
load Reference_data.mat    %GPS��������

disp('Step2:��ʼ������;');
%% һЩ��������������
WIE = 7.292115e-5;           % ������ת���ٶ�
r0  = 6378137.0;             % ����뾶
EE  =  0.0818191908426;      % ƫ����
d2r   =  pi/180;             % degree to radian�Ƕ�ת����
r2d   =  180/pi;             % radian to degree����ת�Ƕ�
dh2rs =  d2r/3600;           % deg/h to rad/s���ٶ�ת���ٶ�
%% ��������ϵ�³�ʼ����̬���ٶȣ�λ��
yaw = (0)*pi/180;%�����
pitch = 0*pi/180;%������
roll = 0*pi/180;%������
cbn=eul2dcm(roll,pitch,yaw);
cnb=cbn';
q=dcm2quat(cbn)';

Vn=0;%�����ٶ�
Ve=0;%�����ٶ�
Vd=0;%�����ٶ�
V_last=[Vn Ve Vd]';

Lati = 31.4913627505302*pi/180;%γ��
Longi= 120.849577188492*pi/180;%����
Alti = 6.6356;%�߶�

sampt0=1/200;%�ߵ�ϵͳ����ʱ��

Rn = r0*(1-EE^2)/(1-EE^2*(sin(Lati))^2)^1.5;         %����Ȧ���ʰ뾶
Re = r0/(1-EE^2*(sin(Lati))^2)^0.5;                  %î��Ȧ���ʰ뾶
g_u = -9.7803267711905*(1+0.00193185138639*sin(Lati)^2)...
    /((1-0.00669437999013*sin(Lati)^2)^0.5 *(1.0 + Alti/r0)^2);
g = [0 0 -g_u]';%����
g0=9.80665;
%% �������˲�P��Q��R����
% P�����ã����Ƶ�Э�������
std_roll    =   (5)*d2r;
std_pitch   =   (5)*d2r;
std_yaw     =   (60)*d2r;

std_vel     =   0.1;
std_pos     =   5;

std_gyro    =   3*0.5*dh2rs;             % �������Ư��0.5��/Сʱ
std_acc     =   3*0.15e-3*g0;            % �ӱ���ƫ0.15mg

Pfilter     =   diag([std_roll^2 std_pitch^2 std_yaw^2 std_vel^2 std_vel^2 std_vel^2 (std_pos/3600/30/57.3)^2 (std_pos/3600/30/57.3)^2 std_pos^2  std_gyro^2 std_gyro^2 std_gyro^2 std_acc^2 std_acc^2 std_acc^2]);

% Q�����ã�����������Э�������
std_Wg      =   0.15*(2.909*1e-4);       % ����Ư������,��/����Сʱת����rad/������
std_Wa      =   0.21/60/3;               % �ӱ�Ư������
Qkf         =   diag([std_Wg^2 std_Wg^2 std_Wg^2 std_Wa^2 std_Wa^2 std_Wa^2]);

G			=   zeros(15, 6);
F = zeros(15);
F_i=zeros(9,9);
F_s=zeros(9,6);
H           =   zeros(6,15);
H(1:3,4:6)  =   eye(3);
H(4:6,7:9)  =   eye(3);
% R�����ã���������Э�������
R           =   diag([std_vel^2 std_vel^2 std_vel^2 (std_pos/3600/30/57.3)^2 (std_pos/3600/30/57.3)^2 (std_pos)^2]);
%��ʼ��һЩ����
cnt_GPS=1;                          % GPS������
k_StartGPS=15;                      % GPS��ʼ����λ�õ�֡��
X_kalman=zeros(size(Pfilter,1),1);  % ״̬��
P_kalman=Pfilter;                   % ����״̬Э�������
output_data=cell(1,length(IMU));    % ����EKF������ĵ�����Ϣ

% ��ʼ����
for cnt_IMU=1:length(IMU) %ѭ�����г���
    if mod(cnt_IMU,5000)==0 %���5000��������н���
       disp(['Time: ', num2str(IMU{cnt_IMU}.time)]);
    end
    
    % ��ȡIMU����
    WM=IMU{cnt_IMU}.gyro'-X_kalman(10:12);%ȥ����ƫ
    AM=IMU{cnt_IMU}.acce'-X_kalman(13:15);%������ʵ����ֵ
    
    % ����ѧ״̬���̵�״̬ת�ƾ���F��������
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
    % ����ѧ״̬ת�ƾ����Э�������Q
    Q_cont      =   Qkf*sampt0;
    Q_kalman =   [[Q_cont, zeros(6,3)];[zeros(3,6), zeros(3,3)]];
    
    % ����״̬��
    X_kalman=F*X_kalman+[zeros(9,1);g];
    % ������̬��update����ʹ��
    [roll, pitch, yaw]=dcm2angle(cbn,'zyx');
    % �;�����̬����
    % ���ٶȻ��ֵõ���̬
    [cbn, q]=gyro_update(cbn, WM, sampt0, q);
    % ���ٶ������ĵ�ͨ�˲��õ���̬
    [cbn, q]=accel_update(cbn, AM, sampt0, q, roll, pitch);
    % ���㿨�����˲���״̬���Ƽ��������    
    P_kalman=F*P_kalman*F'+Q_kalman; %״̬Ԥ��Э�������
    K=P_kalman*H'*inv(H*P_kalman*H'+R);%���㿨�����������
    % ����״̬����
    if cnt_GPS<=length(Reference_data) && IMU{cnt_IMU}.time==Reference_data{cnt_GPS}.time
        % �����ǰʱ����GPS�۲�ֵ
        [vn,ve,vr,LL]=BLH2xyzV(Reference_data{cnt_GPS}.latitude*pi/180,Reference_data{cnt_GPS}.longitude*pi/180,Reference_data{cnt_GPS}.altitude,0,0,0,0,0,0,0,0,jd);%�����ο�λ��
        z_gps=[vn-IMU{cnt_IMU}.vn; ve-IMU{cnt_IMU}.ve; vr-IMU{cnt_IMU}.vd;LL(1)*r2d*3600;LL(2)*r2d*3600;LL(3)];
        X_kalman=X_kalman+K*(z_gps(1:6)-H(1:6,:)*X_kalman);%״̬��������
        P_kalman=(eye(size(K,1))-K*H)*P_kalman;%Э����������
        cnt_GPS=cnt_GPS+1;%GPS����������
    else
        % �����ǰʱ��û��GPS�۲�ֵ
        X_kalman=X_kalman;%״̬��������
        P_kalman=P_kalman;%Э����������
    end
    
    %�������е�����Ϣ
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