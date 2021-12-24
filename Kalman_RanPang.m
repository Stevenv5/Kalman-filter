clc 
clear all
close  all;

%% 迭代设定
delta_t=1;%采样周期   
step=1:delta_t:30;%迭代第一步，采样周期和最后一步
N =length(step);%迭代次数

%%  矩阵的设定
A=[1 1;0 1];%状态转移矩阵（表示状态变量从上一步到这一步的变化）
B=[]; %控制输入矩阵
H=[1 0;0 1];%实际测量的状态转移矩阵（表示状态变量到测量变量的转移过程）

Q=[0.1 0;0 0.1];%过程噪声（表示模型的不确定性）
R=[1 0;0 1];%测量噪声（表示测量的不确定性）

% 初始化
X0=[0;1];%初始位置和速度
P0=[1,0;0,1];%状态估计误差的协方差矩阵初始化
X_true=X0; %状态实际值初始化
X_posterior=X0;%后验估计初始化
P_posterior=P0;%后验状态估计误差的协方差矩阵初始化

%开始记录状态变量
speed_true = [];
position_true = [];

speed_measure = [];
position_measure = [];

speed_prior_est = [];
position_prior_est = [];

speed_posterior_est = [];
position_posterior_est = [];

%%Kalman filter的预测和校正的迭代过程
for i=1:N
    %生成真实值并保存
    w=sqrt(0.1)*randn(2,1);%生成方差（统计学中的方差相当于工程上的误差）为0.1,均值为0，服从高斯分布的过程噪声（2*1维）
    X_true=A*X_true+w;
    position_true(i,1)=X_true(1); %保存这一步的真实位置
    speed_true(i,1)=X_true(2); %保存这一步的真实速度
    
    %生成测量值
    v=sqrt(1)*randn(2,1); %生成方差为1,均值为0，服从高斯分布的测量噪声（2*1维）
    Z_measure=H*X_true+v; %状态变量的测量值=真值+测量噪声（这是仿真数据，实际中应直接通过实验采集的数据值作为测量值）
    position_measure(i,1)=X_true(1);
    speed_measure(i,1)=X_true(2);
    
    %先验估计
    X_prior=A*X_posterior; %①
    position_prior_est(i,1)=X_prior(1);
    speed_prior_est(i,1)=X_prior(2);
    
    %计算先验状态估计误差的协方差矩阵
    P_prior=A*P_posterior*A'+Q ; %②
    
    %计算卡尔曼增益
    K = (P_prior * H')/(H * P_prior * H' + R);  %③
    
    %后验估计
    X_posterior=X_prior+K*(Z_measure-H*X_prior);  %④
    position_posterior_est(i,1)=X_posterior(1);
    speed_posterior_est(i,1)=X_posterior(2);
    
    %更新状态估计误差的协方差矩阵
    Pk = (eye(2)-K*H)*P_prior;  %⑤
end




%% 可视化
figure %位置可视化
plot(step,position_true,'-^k','linewidth',1);
hold on
plot(step,position_measure,'--or','linewidth',1);
hold on
plot(step,position_prior_est,':+g','linewidth',1);
hold on
plot(step,position_posterior_est,'-.*b','linewidth',1);
legend('true','measure','prior est','posterior est','Location','Northwest');
xlabel('t(s)')
ylabel('S(m)')
title('Prediction position of Kalman filter');

figure%速度可视化
plot(step,speed_true,'-^k','linewidth',1);
hold on
plot(step,speed_measure,'--or','linewidth',1);
hold on
plot(step,speed_prior_est,':+g','linewidth',1);
hold on
plot(step,speed_posterior_est,'-.*b','linewidth',1);
legend('true','measure','prior est','posterior est','Location','Southwest');
xlabel('t(s)')
ylabel('V(m/s)')
title('Prediction speed of Kalman filter');