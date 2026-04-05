function Power_adjust = power_NREM_nor(NREM_duration, power_NE)
%ESTIMATE_POWER_YOUNG 计算基于调整后的线性模型的 Power_Young
%
% 输入:
%   Young_NREM_baseline - NREM duration vector（必须为列向量）
%   Young_power_baseline - Power vector（必须为列向量）
%
% 输出:
%   Power_Young - 调整后估计的 power（已从 ln 空间还原）

    % 去除非法值（<=0 或 NaN）
    validIdx = NREM_duration > 0 & power_NE > 0 & ...
               ~isnan(NREM_duration) & ~isnan(power_NE);

    x = log(NREM_duration(validIdx));
    y = log(power_NE(validIdx));

    % 中心化处理
    x_adj = x - 5.81;
    y_adj = y - 1.78;

    % 计算每个点的 slope a
    a = y_adj ./ x_adj;

    % 线性拟合反算，得到 ln(power)
    Power_adjusted = a .* 1 + (1.78 - 5.81 .* a);

    % 指数还原
    Power_adjust = exp(Power_adjusted);
end