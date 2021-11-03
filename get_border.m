function [leftpoint, rightpoint] = get_border(x0, y0, x1, y1, w_left0, w_left1, w_right0, w_right1, oldleft, oldright)
    A = [x0 y0]; 
    B = [x1 y1]; 
    Clen = (w_left0+w_left1 + w_right0+w_right1)/2;
    slope = (B(2)-A(2)) / (B(1)-A(1)); 
    yint = B(2) - slope*B(1);
    C(1) = range([B(1),A(1)])/2+min([A(1),B(1)]); 
    C(2) = slope * C*(1) + yint;
    perpSlope = -1/slope; 
    %perpYint = C(2) - perpSlope*C(1);
    x = C(1) + (Clen*sqrt(1/(1+perpSlope^2)))*[-1,1]; 
    y = C(2) + (perpSlope*Clen*sqrt(1/(1+perpSlope^2)))*[-1,1];
    if oldleft(1) == 0 && oldright(1) == 0
        leftpoint = [x(1), y(1)];
        rightpoint = [x(2), y(2)];
    else
        disp((x(1)-oldleft(1))^2 + (y(1)-oldleft(2))^2)
        disp( (x(1)-oldright(1))^2 + (y(1)-oldright(2))^2)
        if (x(1)-oldleft(1))^2 + (y(1)-oldleft(2))^2 < (x(1)-oldright(1))^2 + (y(1)-oldright(2))^2
            disp("Left")
            leftpoint = [x(1), y(1)];
            rightpoint = [x(2), y(2)];
        else
            disp("Right")
            leftpoint = [x(2), y(2)];
            rightpoint = [x(1), y(1)];
        end
    end
end