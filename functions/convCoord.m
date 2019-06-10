function [x,y] = convCoord(c_x,c_y,winFE,UFE)

lengh_tot_x=UFE(1,end)-UFE(1,1);
lengh_x=UFE(1,end)-c_x;

lengh_tot_y=UFE(2,end)-UFE(2,1);
lengh_y=UFE(2,end)-c_y;

sizeFE = abs(winFE(1,:)-winFE(2,:))+1;

x=sizeFE(1)-(lengh_x*(sizeFE(1)-1))/lengh_tot_x;
% x=round(x);
x=round(x)+1;

y=sizeFE(2)-(lengh_y*(sizeFE(2)-1))/lengh_tot_y;
% y=round(y);
y=round(y)+1;