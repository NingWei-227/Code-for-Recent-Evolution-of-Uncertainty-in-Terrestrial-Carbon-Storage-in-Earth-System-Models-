function tuaE_scalar=cal_scaler(x,inputs)
    Q10=x(1);
    Tem=inputs(1,:);
    maxTem=max(Tem);
    pre=inputs(2,:);
    maxPre =max(pre);
    scaler_tem=Q10.^((Tem-maxTem)/10);
    scaler_pre=pre/maxPre;
    T_scaler=scaler_tem.*scaler_pre;
    tuaE_scalar=[scaler_tem;scaler_pre;T_scaler];
end