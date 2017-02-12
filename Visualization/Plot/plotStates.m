
% plot states

hold on
if plot_flag(1)
    plot(T,P(:,1),'b','DisplayName','x (m)') 
end
if plot_flag(2)
    plot(T,P(:,2),'r','DisplayName','y (m)')
end
if plot_flag(3)
    plot(T,P(:,3),'Color',[0 0.5 0],'DisplayName','phi (rad)')
end
if plot_flag(4)
    plot(T,P(:,4),'k','DisplayName','alpha (rad)')
end
if plot_flag(5)
    plot(T,P(:,5),'m','DisplayName','beta (rad)')
end
if plot_flag(6)
    plot(T,P(:,6),'b--','DisplayName','xDot (m/s)')
end
if plot_flag(7)
    plot(T,P(:,7),'r--','DisplayName','yDot (m/s)')
end
if plot_flag(8)
    plot(T,P(:,8),'--','Color',[0 0.5 0],'DisplayName','phiDot (rad/s)')
end
if plot_flag(9)
    plot(T,P(:,9),'k--','DisplayName','alphaDot (rad/s)')
end
if plot_flag(10)
    plot(T,P(:,10),'m--','DisplayName','betaDot (rad/s)')
end
if plot_flag(11)
    plot(T,P(:,11),'Color',[0.5 0 0],'DisplayName','L (m)')
    plot([T(1) T(end)],[Yparam.L_sp0 Yparam.L_sp0],'k--','DisplayName','L_{eq,R} (m)')
end
if plot_flag(12)
    plot(T,P(:,12),'--','Color',[0.5 0 0],'DisplayName','LDot_R (m/s)')
end
if plot_flag(13)
    plot(T,P(:,13),'--','Color',[0.8 0 0],'DisplayName','Energy (J)')
end
if plot_flag(14)
    plot(T,P(:,14),'--','Color',[0 0 0.8],'DisplayName','Desired Energy (J)')
end
if plot_flag(15)
    plot(T,P(:,15),'--','Color',[0.8 0 0],'DisplayName','tau_{hip,R} (N*m)')
end
if plot_flag(16)
    plot(T,P(:,16),'--','Color',[0 0 0.8],'DisplayName','tau_{knee,R} (N*m)')
end
if plot_flag(17)
    plot(T,P(:,17),'Color',[1 0 0],'DisplayName','F_cx (N)')
end
if plot_flag(18)
    plot(T,P(:,18),'Color',[0 0 1],'DisplayName','F_cy (N)')
end
if plot_flag(19)
    plot(T,P(:,19),'Color',[1 0 0],'DisplayName','theta (rad)')
end
if plot_flag(20)
    plot(T,P(:,20),'Color',[0 0 1],'DisplayName','dtheta (rad/s)')
end
if plot_flag(21)
    plot(T,P(:,21),'Color',[0 0 1],'DisplayName','RightFootPos_x (m)')
end
if plot_flag(22)
    plot(T,P(:,22),'Color',[0 0 1],'DisplayName','RightFootPos_y (m)')
end
if plot_flag(23)
    plot(T,P(:,23),'--','Color',[0.8 0 0],'DisplayName','tau_{hip,L} (N*m)')
end
if plot_flag(24)
    plot(T,P(:,24),'--','Color',[0 0 0.8],'DisplayName','tau_{knee,L} (N*m)')
end
hold off
