% 用四阶离散格式去近似积分结果;

% 这里的f_matrix是一个空间剖分行，时间剖分列的矩阵;
% 时间步长dt;
function integral = integral_appro(f_matrix, left, right, dt)

    % 考虑零这个点的存在，matlab从一开始编号;
    left = floor(left/dt) + 1;
    right = floor(right/dt) + 1;

    % 计算求和项;
    f_sum = f_matrix(:,:,left);
    for i = left+1:right
        f_sum = f_sum + f_matrix(:,:,i);
    end

    % 计算积分;
    integral = -5/8*(f_matrix(:,:,left) + f_matrix(:,:,right));
    integral = integral + 1/6*(f_matrix(:,:,left+1) + f_matrix(:,:,right-1));
    integral = integral - 1/24*(f_matrix(:,:,left+2) + f_matrix(:,:,right-2));
    integral = integral + f_sum;
    integral = dt*integral;
    
end
    
