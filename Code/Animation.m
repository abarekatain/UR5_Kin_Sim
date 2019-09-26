
% jnt1 = [0; 0; 89.2];
% jnt2 = [-425.0.*cos(q1(i)).*cos(q2(i)); -425.0.*cos(q2(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 89.2];
% jnt3 = [392.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i)) - 392.0.*cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 425.0.*cos(q1(i)).*cos(q2(i)); 392.0.*sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 425.0.*cos(q2(i)).*sin(q1(i)) - 392.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 392.0.*cos(q2(i)).*sin(q3(i)) + 392.0.*cos(q3(i)).*sin(q2(i)) + 89.2];
% jnt4 = [392.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i)) - 425.0.*cos(q1(i)).*cos(q2(i)) - 392.0.*cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 109.3.*sin(q1(i)); 109.3.*cos(q1(i)) - 425.0.*cos(q2(i)).*sin(q1(i)) + 392.0.*sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 392.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 392.0.*cos(q2(i)).*sin(q3(i)) + 392.0.*cos(q3(i)).*sin(q2(i)) + 89.2];
% jnt5 = [94.75.*sin(q4(i)).*(cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 1.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i))) - 425.0.*cos(q1(i)).*cos(q2(i)) - 109.3.*sin(q1(i)) + 94.75.*cos(q4(i)).*(cos(q1(i)).*cos(q2(i)).*sin(q3(i)) + cos(q1(i)).*cos(q3(i)).*sin(q2(i))) - 392.0.*cos(q1(i)).*cos(q2(i)).*cos(q3(i)) + 392.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i)); 109.3.*cos(q1(i)) - 425.0.*cos(q2(i)).*sin(q1(i)) - 94.75.*sin(q4(i)).*(sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 1.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i))) + 94.75.*cos(q4(i)).*(cos(q2(i)).*sin(q1(i)).*sin(q3(i)) + cos(q3(i)).*sin(q1(i)).*sin(q2(i))) + 392.0.*sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 392.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 392.0.*cos(q2(i)).*sin(q3(i)) + 392.0.*cos(q3(i)).*sin(q2(i)) + 94.75.*cos(q4(i)).*(cos(q2(i)).*cos(q3(i)) - 1.0.*sin(q2(i)).*sin(q3(i))) - 94.75.*sin(q4(i)).*(cos(q2(i)).*sin(q3(i)) + cos(q3(i)).*sin(q2(i))) + 89.2];
% jnt6 = [94.75.*sin(q4(i)).*(cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 1.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i))) - 425.0.*cos(q1(i)).*cos(q2(i)) - 82.5.*sin(q5(i)).*(cos(q4(i)).*(cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 1.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i))) - 1.0.*sin(q4(i)).*(cos(q1(i)).*cos(q2(i)).*sin(q3(i)) + cos(q1(i)).*cos(q3(i)).*sin(q2(i)))) - 82.5.*cos(q5(i)).*sin(q1(i)) - 109.3.*sin(q1(i)) + 94.75.*cos(q4(i)).*(cos(q1(i)).*cos(q2(i)).*sin(q3(i)) + cos(q1(i)).*cos(q3(i)).*sin(q2(i))) - 392.0.*cos(q1(i)).*cos(q2(i)).*cos(q3(i)) + 392.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i)); 109.3.*cos(q1(i)) + 82.5.*cos(q1(i)).*cos(q5(i)) - 425.0.*cos(q2(i)).*sin(q1(i)) - 94.75.*sin(q4(i)).*(sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 1.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i))) + 82.5.*sin(q5(i)).*(cos(q4(i)).*(sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 1.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i))) + sin(q4(i)).*(cos(q2(i)).*sin(q1(i)).*sin(q3(i)) + cos(q3(i)).*sin(q1(i)).*sin(q2(i)))) + 94.75.*cos(q4(i)).*(cos(q2(i)).*sin(q1(i)).*sin(q3(i)) + cos(q3(i)).*sin(q1(i)).*sin(q2(i))) + 392.0.*sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 392.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 392.0.*cos(q2(i)).*sin(q3(i)) + 392.0.*cos(q3(i)).*sin(q2(i)) + 94.75.*cos(q4(i)).*(cos(q2(i)).*cos(q3(i)) - 1.0.*sin(q2(i)).*sin(q3(i))) - 94.75.*sin(q4(i)).*(cos(q2(i)).*sin(q3(i)) + cos(q3(i)).*sin(q2(i))) + 82.5.*sin(q5(i)).*(sin(q4(i)).*(cos(q2(i)).*cos(q3(i)) - 1.0.*sin(q2(i)).*sin(q3(i))) + cos(q4(i)).*(cos(q2(i)).*sin(q3(i)) + cos(q3(i)).*sin(q2(i)))) + 89.2];

base = [0;0;0];

t = tout;
animate_step = 20;

for i = 1:animate_step:size(t)
    
%    if i==1 
%        pause; 
%    end
   
   clf
   hold;
   
   axis([-500 500 -500 500 -500 500 -500 500])
   
%    plot3(x_des(:,1),x_des(:,2),x_des(:,3),'Color','black');
   
   jnt1 = [0; 0; 89.2];
   jnt2 = [-425.0.*cos(q1(i)).*cos(q2(i)); -425.0.*cos(q2(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 89.2];
   jnt3 = [392.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i)) - 392.0.*cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 425.0.*cos(q1(i)).*cos(q2(i)); 392.0.*sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 425.0.*cos(q2(i)).*sin(q1(i)) - 392.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 392.0.*cos(q2(i)).*sin(q3(i)) + 392.0.*cos(q3(i)).*sin(q2(i)) + 89.2];
   jnt4 = [392.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i)) - 425.0.*cos(q1(i)).*cos(q2(i)) - 392.0.*cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 109.3.*sin(q1(i)); 109.3.*cos(q1(i)) - 425.0.*cos(q2(i)).*sin(q1(i)) + 392.0.*sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 392.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 392.0.*cos(q2(i)).*sin(q3(i)) + 392.0.*cos(q3(i)).*sin(q2(i)) + 89.2];
   jnt5 = [94.75.*sin(q4(i)).*(cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 1.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i))) - 425.0.*cos(q1(i)).*cos(q2(i)) - 109.3.*sin(q1(i)) + 94.75.*cos(q4(i)).*(cos(q1(i)).*cos(q2(i)).*sin(q3(i)) + cos(q1(i)).*cos(q3(i)).*sin(q2(i))) - 392.0.*cos(q1(i)).*cos(q2(i)).*cos(q3(i)) + 392.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i)); 109.3.*cos(q1(i)) - 425.0.*cos(q2(i)).*sin(q1(i)) - 94.75.*sin(q4(i)).*(sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 1.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i))) + 94.75.*cos(q4(i)).*(cos(q2(i)).*sin(q1(i)).*sin(q3(i)) + cos(q3(i)).*sin(q1(i)).*sin(q2(i))) + 392.0.*sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 392.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 392.0.*cos(q2(i)).*sin(q3(i)) + 392.0.*cos(q3(i)).*sin(q2(i)) + 94.75.*cos(q4(i)).*(cos(q2(i)).*cos(q3(i)) - 1.0.*sin(q2(i)).*sin(q3(i))) - 94.75.*sin(q4(i)).*(cos(q2(i)).*sin(q3(i)) + cos(q3(i)).*sin(q2(i))) + 89.2];
   jnt6 = [94.75.*sin(q4(i)).*(cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 1.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i))) - 425.0.*cos(q1(i)).*cos(q2(i)) - 82.5.*sin(q5(i)).*(cos(q4(i)).*(cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 1.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i))) - 1.0.*sin(q4(i)).*(cos(q1(i)).*cos(q2(i)).*sin(q3(i)) + cos(q1(i)).*cos(q3(i)).*sin(q2(i)))) - 82.5.*cos(q5(i)).*sin(q1(i)) - 109.3.*sin(q1(i)) + 94.75.*cos(q4(i)).*(cos(q1(i)).*cos(q2(i)).*sin(q3(i)) + cos(q1(i)).*cos(q3(i)).*sin(q2(i))) - 392.0.*cos(q1(i)).*cos(q2(i)).*cos(q3(i)) + 392.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i)); 109.3.*cos(q1(i)) + 82.5.*cos(q1(i)).*cos(q5(i)) - 425.0.*cos(q2(i)).*sin(q1(i)) - 94.75.*sin(q4(i)).*(sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 1.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i))) + 82.5.*sin(q5(i)).*(cos(q4(i)).*(sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 1.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i))) + sin(q4(i)).*(cos(q2(i)).*sin(q1(i)).*sin(q3(i)) + cos(q3(i)).*sin(q1(i)).*sin(q2(i)))) + 94.75.*cos(q4(i)).*(cos(q2(i)).*sin(q1(i)).*sin(q3(i)) + cos(q3(i)).*sin(q1(i)).*sin(q2(i))) + 392.0.*sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 392.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 392.0.*cos(q2(i)).*sin(q3(i)) + 392.0.*cos(q3(i)).*sin(q2(i)) + 94.75.*cos(q4(i)).*(cos(q2(i)).*cos(q3(i)) - 1.0.*sin(q2(i)).*sin(q3(i))) - 94.75.*sin(q4(i)).*(cos(q2(i)).*sin(q3(i)) + cos(q3(i)).*sin(q2(i))) + 82.5.*sin(q5(i)).*(sin(q4(i)).*(cos(q2(i)).*cos(q3(i)) - 1.0.*sin(q2(i)).*sin(q3(i))) + cos(q4(i)).*(cos(q2(i)).*sin(q3(i)) + cos(q3(i)).*sin(q2(i)))) + 89.2];

   
   plot3([base(1) jnt1(1)], [base(2) jnt1(2)], [base(3) jnt1(3)],'Marker','o','LineWidth',3);
   plot3([jnt1(1) jnt2(1)], [jnt1(2) jnt2(2)], [jnt1(3) jnt2(3)],'Marker','o','LineWidth',3);
   plot3([jnt2(1) jnt3(1)], [jnt2(2) jnt3(2)], [jnt2(3) jnt3(3)],'Marker','o','LineWidth',3);
   plot3([jnt3(1) jnt4(1)], [jnt3(2) jnt4(2)], [jnt3(3) jnt4(3)],'Marker','o','LineWidth',3);
   plot3([jnt4(1) jnt5(1)], [jnt4(2) jnt5(2)], [jnt4(3) jnt5(3)],'Marker','o','LineWidth',3);
   plot3([jnt5(1) jnt6(1)], [jnt5(2) jnt6(2)], [jnt5(3) jnt6(3)],'Marker','o','LineWidth',3);
   
   
   pause(0.1)
end


for i = 1:animate_step:size(t)
    hold on
    x = [(379*cos(q2(i) + q3(i))*cos(q1(i))*sin(q4(i)))/4 - 425*cos(q1(i))*cos(q2(i)) - (165*cos(q5(i))*sin(q1(i)))/2 - (165*cos(q2(i) + q3(i) + q4(i))*cos(q1(i))*sin(q5(i)))/2 - (1093*sin(q1(i)))/10 + (379*sin(q2(i) + q3(i))*cos(q1(i))*cos(q4(i)))/4 - 392*cos(q1(i))*cos(q2(i))*cos(q3(i)) + 392*cos(q1(i))*sin(q2(i))*sin(q3(i));
        (1093*cos(q1(i)))/10 + (165*cos(q1(i))*cos(q5(i)))/2 - 425*cos(q2(i))*sin(q1(i)) + 392*sin(q1(i))*sin(q2(i))*sin(q3(i)) - (165*cos(q2(i) + q3(i) + q4(i))*sin(q1(i))*sin(q5(i)))/2 + (379*cos(q2(i) + q3(i))*sin(q1(i))*sin(q4(i)))/4 + (379*sin(q2(i) + q3(i))*cos(q4(i))*sin(q1(i)))/4 - 392*cos(q2(i))*cos(q3(i))*sin(q1(i));
        392*sin(q2(i) + q3(i)) + 425*sin(q2(i)) + sin(q5(i))*((165*cos(q2(i) + q3(i))*sin(q4(i)))/2 + (165*sin(q2(i) + q3(i))*cos(q4(i)))/2) + (379*cos(q2(i) + q3(i))*cos(q4(i)))/4 - (379*sin(q2(i) + q3(i))*sin(q4(i)))/4 + 446/5;
        0;0;0];
    plot3(x(1),x(2),x(3),'r--','Color','red')
end
