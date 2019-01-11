% 自洽场柔性链问题;
% 二维问题;
clear; clc; close all;
delete('./result/*.txt');
delete('./figure/*.eps');
delete('./figure/*.png');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%         参数           %%%%%%%%%%%%%%%%%%%%%%

% 参数设置;
chiN = 25;
f = 0.2;
L = 20;

% 柔性链的数目;
n = 1;

% 设置空间剖分;
N = 128;
r=linspace(0,L,N);
% 时间步长dt;
dt = 0.01                                                                                                                ;
dt1 = 2.0;
dt2 = 2.0;

% 设置迭代终止的误差范数;
norm_eps = 10^-6;

% 设置空间维度;
ncpt = N*ones(1, 2);

% 设置复空间中的k^2值;
ksquare = zeros(ncpt);
for j1 = 1:ncpt(1)
    if (j1 > N/2+1) k1 = j1-N;
    else k1 = j1;
    end
    for j2 = 1:ncpt(2)
        if (j2 > N/2+1) k2 = j2-N;
        else k2 = j2;
        end
        k = [k1, k2];
        k = k-1;
        % 计算k^2;
        ksquare(j1,j2) = sum(k.^2);
    end
end
ksquare = (2*pi/L)^2*ksquare;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%         初值           %%%%%%%%%%%%%%%%%%%%%%

% 实空间的mu+初值;
mu_plus = zeros(ncpt);
% 得到mu-的初值;
% mu_fft_minus = zeros(ncpt);
mu_minus = zeros(ncpt);
% %层状相
% for i = 1:ncpt(1)
%     mu_minus(:,i)=sin(4*r(i));
% end
%六状相
for i = 1:ncpt(1)
    for j = 1:ncpt(2)
       mu_minus(j,i) = cos(cos(cos(pi/3*0)*r(i)+sin(pi/3*0)*r(j)))+cos(cos(cos(pi/3*1)*r(i)+sin(pi/3*1)*r(j)))+cos(cos(cos(pi/3*2)*r(i)+sin(pi/3*2)*r(j)))+cos(cos(cos(pi/3*3)*r(i)+sin(pi/3*3)*r(j)))+cos(cos(cos(pi/3*4)*r(i)+sin(pi/3*4)*r(j)))+cos(cos(cos(pi/3*5)*r(i)+sin(pi/3*5)*r(j)));
    end
end
% 得到对应索引;
%kindex = [  0,   1;
%             0,  -1;
%             1,   0;
%            -1,   0;
%             1,   1;
%             1,  -1;
%            -1,   1;
%            -1,  -1];
% 
% [nr, nc] = size(kindex);

%%%%%%%%%%%% 将复空间索引变成matlab编号;
% for i = 1:nr
%     for j = 1:nc
%         if (kindex(i,j) < 0) 
%             kindex(i,j) = kindex(i,j)+N;
%         end
%     end
% end
% kindex = kindex+1;
% 
% % 赋非零值;
% for i = 1:nr
%     mu_fft_minus(kindex(i,1),kindex(i,2)) = 0.1;
% end
% 
% % 变换到实空间;
% mu_minus = fftn(mu_fft_minus);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%         迭代           %%%%%%%%%%%%%%%%%%%%%%

% 输出重定向;
diary('./result/flexible.txt');
diary on;


% 初始化误差范数;
% mu_plus_norm = 1;
% mu_minus_norm = 1;
ediff=inf;
res=inf;
% 初始化计数器;
num = 0;

% 存储mu+和mu-的向量;
mu_plus_norm_vector = [0];
mu_minus_norm_vector = [0];

% 存储Hamilt的向量;
Hamilt_vector = [0];
Hold = inf;
% 进行循环迭代;
while (res > norm_eps)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%        解PDE           %%%%%%%%%%%%%%%%%%%%%%
    res=abs(ediff);
    
    % 计算omega;
    omega_a = mu_plus-mu_minus;
    omega_b = mu_plus+mu_minus;

    % 解关于q的PDE系统;
    q_left = PDE_solve(ones(ncpt), omega_a, 0, f, ksquare, dt);
    q_right = PDE_solve(q_left(:,:,end), omega_b, f, 1, ksquare, dt);

    % 组装q;
    q_right = q_right(:,:,2:1:end);
    q_matrix = cat(3, q_left, q_right);
    

    % 解关于q_plus的PDE系统;
    q_plus_right = PDE_solve(ones(ncpt), omega_b, 0, 1-f, ksquare, dt);
    q_plus_left = PDE_solve(q_plus_right(:,:,end), omega_a, 1-f, 1, ksquare, dt);

    % 组装q_plus;
    % 矩阵倒装;
    q_plus_right = q_plus_right(:,:,end:-1:1);
    q_plus_left = q_plus_left(:,:,end:-1:2);
    q_plus_matrix = cat(3, q_plus_left, q_plus_right);
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%        Q,phi           %%%%%%%%%%%%%%%%%%%%%%

    % 计算Q;
    % 取t=1时刻的第零个Fourier系数;
    Q = ifftn(q_matrix(:,:,end));
    Q = real(Q(1));
    fprintf('%15f',Q)
    % 计算phi_a, phi_b;
    q_mul_plus = q_matrix.*q_plus_matrix;

    % 计算phi_a, phi_b;
    phi_a = 1/Q*integral_appro(q_mul_plus, 0, f, dt);
    phi_b = 1/Q*integral_appro(q_mul_plus, f, 1, dt);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%  更新mu_plus, mu_minus %%%%%%%%%%%%%%%%%%%%%%

    % 一阶时间离散计算mu_plus;
    mu_plus = mu_plus + (phi_a+phi_b-1)*dt1;

    % 一阶时间离散计算mu_minus;
    mu_minus = mu_minus + (phi_a-phi_b)*dt2;
    mu_minus = 1/(1+2*dt2/chiN)*mu_minus;

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%       计算误差范数     %%%%%%%%%%%%%%%%%%%%%%

    % 计算H关于mu_plus和mu_minus的变份;
    mu_plus_norm = max(max(abs(phi_a+phi_b-1)));
    mu_minus_norm = max(max(abs(phi_a-phi_b-2/chiN*mu_minus)));

    % 存储误差范数;
    %mu_plus_norm_vector(end+1) = mu_plus_norm;
    %mu_minus_norm_vector(end+1) = mu_minus_norm;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%       计算能量泛函     %%%%%%%%%%%%%%%%%%%%%%

    % Hamilt的第一个部分,\int d\bm{r} (-\mu_{+} + \frac{1}{\chi N} \mu_{-}^2);
    Hamilt1 = -mu_plus + 1/chiN*mu_minus.^2;
    Hamilt1 = ifftn(Hamilt1);
    Hamilt1 = real(Hamilt1(1));

    % Hamilt的第二个部分;
    Hamilt2 = -log(Q);

    % 总的Hamilt能量;
    Hamilt = Hamilt1 + Hamilt2;
    Hamilt = n*Hamilt;
    fprintf('%15f',Hamilt);
    % 存储总能量;
    Hamilt_vector(end+1) = Hamilt;
    
    ediff = Hamilt - Hold;
    Hold = Hamilt;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%        格式化输出      %%%%%%%%%%%%%%%%%%%%%%

    % 计数器迭代;
    num = num + 1;
    
    % 输出;
    fprintf('迭代步长：%d  误差范数：%.15e \t Hamilt：%.15e\n\n',...
        num, ediff, Hamilt);
    
%     fprintf('迭代步长：%d  误差范数：+：%.6e  -：%.6e \t Hamilt：%.12e\n\n',...
%         num, mu_plus_norm, mu_minus_norm, Hamilt);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%       检验Q值           %%%%%%%%%%%%%%%%%%%%%%


length_mul = length(q_mul_plus(1,1,:));

for i = 1:length_mul
    Q = ifftn(q_mul_plus);
    Q = Q(1);
    fprintf('时间s=%.6f \t 能量Q=%.6f\n',(i-1)*dt,Q);
end
imagesc(phi_a);
colorbar;
% axis off;

% 结束输出重定向;
diary off;


% 保存数据;
save('./result/flexible_phi.mat', 'phi_a', 'phi_b');
save('./result/flexible_error.mat', 'mu_plus_norm_vector', 'mu_plus_norm_vector');
save('./result/flexible_hamilt.mat', 'Hamilt_vector');

