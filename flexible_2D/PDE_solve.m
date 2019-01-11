% 解PDE系统，得到q_plus的嵌段;
% 算子分裂;

% mu_plus表示mu+; mu_minus表示mu-;
% left是迭代起始时刻; right是迭代终止时刻;
% ksquare是Fourier空间的k^2;
% dt是时间步长;

function q_matrix = PDE_solve(q_initial, omega, left, right, ksquare, dt)

    % 计算对应的系数;
    % 应用全隐格式;
    % 第一步和第三步的系数(半步);
    coeff_fftn = exp(-dt/2*ksquare);
    % 第二步的系数(一步);
    coeff_real = exp(-dt*omega); 

    % 初始化q_matrix;
    q_matrix(:,:,1) = q_initial;
    
    % 赋循环初值;
    q_fft_mid = ifftn(q_initial);

    % 进行循环迭代;
    for t = left+dt:dt:right
        % 第一步;
        q_fft_mid  = coeff_fftn.*q_fft_mid;
        q_mid = fftn(q_fft_mid);

        % 第二步;
        q_mid = coeff_real.*q_mid;
        
        % 第三步;
        q_fft_mid = ifftn(q_mid);
        q_fft_mid = coeff_fftn.*q_fft_mid;
        q_matrix(:,:,end+1) = real(fftn(q_fft_mid));
    end

end
