%thrust schedule is a 2×n matrix
function mdot = emdot(t,thrustschedule,m_roc,m_burnout,meco)
    if(m_roc>m_burnout && t<thrustschedule(1,size(thrustschedule,2)) && meco==0)
        mdot=interp1(thrustschedule(1,:),thrustschedule(2,:),t);
    else
        mdot=0;
    end
        
