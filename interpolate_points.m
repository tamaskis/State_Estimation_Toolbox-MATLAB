function [x_out,y_out] = interpolate_points(point_1,point_2,p)
    
    x1 = point_1(1);
    y1 = point_1(2);

    x2 = point_2(1);
    y2 = point_2(2);

    theta = atan2(y2-y1,x2-x1);
    
    l = sqrt((x2-x1)^2+(y2-y1)^2);
    
    x_out = x1+p*l*cos(theta);
    y_out = y1+p*l*sin(theta);

end