function fun_min=rmseR2_tuaE(x,inputs)
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
    
    r2=1-sum((tuaE-baseTuaE./T_scaler).^2)./(sum((tuaE-baseTuaE).^2));
    rmse=((sum((tuaE-baseTuaE./T_scaler).^2))./num).^(0.5);
    
    fun_min=rmse/r2;
end