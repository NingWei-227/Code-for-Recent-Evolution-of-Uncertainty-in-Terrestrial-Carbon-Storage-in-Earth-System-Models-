function [cons,ceq]=myfunc_tuaE(x,inputs)
% annual temperature: inputs(1,:)
% annual precipitation: inputs(2,:)
% the modeled tuaE: inputs(3,:)
% set range for R2: 0< R2 <1;

    Q10=x(1);
    baseTuaE=x(2);
    Tem=inputs(1,:);
    maxTem=max(Tem);
    pre=inputs(2,:);
    maxPre =max(pre);
    tuaE=inputs(3,:);
    scaler_tem=Q10.^((Tem-maxTem)/10);
    scaler_pre=pre/maxPre;
    T_scaler=scaler_tem.*scaler_pre;
    
    r2=1-sum((tuaE-baseTuaE./T_scaler).^2)./(sum((tuaE-baseTuaE).^2));

    cons=[r2-1;-r2]; % the inequalities: R2-1<0 and -R2 <0
    ceq=[]; % set the equality []
end