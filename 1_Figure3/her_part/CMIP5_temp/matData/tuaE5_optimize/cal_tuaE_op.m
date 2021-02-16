function tuaE_op= cal_tuaE_op(x,inputs)
% annual temperature: inputs(1,:)
% annual precipitation: inputs(2,:)
% the modeled tuaE: inputs(3,:)
    Q10=x(1);
    baseTuaE=x(2);
    Tem=inputs(1,:);
    maxTem=max(Tem);
    pre=inputs(2,:);
    maxPre =max(pre);
    num=size(inputs,2);
    tuaE=inputs(3,:);
    scaler_tem=Q10.^((Tem-maxTem)/10);
    scaler_pre=pre/maxPre;
    T_scaler=scaler_tem.*scaler_pre;
%      r2=1-sum((values-basedValue./total_s).^2)./(sum((values-basedValue).^2));
%      rmse=((sum((values-basedValue./total_s).^2))./num).^(0.5);
    tuaE_op = baseTuaE./T_scaler;
end