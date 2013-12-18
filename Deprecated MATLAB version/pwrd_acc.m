% externally: keep track of m and dm_dt - need another rk4

function a = pwrd_acc(x,v,m_obj,m_pri,C_D,A_ref,rho_atm,thrust)
    G=6.67384e-11;
    mu=G*m_pri;
    
    r=norm(x);
    a_g=-(x/r)*(mu/(r*r));
    
    a_d=-(v/norm(v))*(1/m_obj)*((1/2)*rho_atm*norm(v)*norm(v)*C_D*A_ref);
    
    a=a_g+a_d+thrust./m_obj;