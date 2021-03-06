function [left_border, right_border] = read_track(track_path)
% MYCOORDINATES computes the coordinates x and y of the points through which the car has to pass at each step. 
%  INPUTS:  
%           track_path   = path to the file .csv containing the track
%
%  OUTPUTS:
%           left_border             = x,y of left border
%           right_border           = x,y of right border

    track = readtable(track_path);
    n = height(track);
    leftpoints = zeros(n-1, 2);
    rightpoints =  zeros(n-1, 2);
    oldleft = [0,0];
    oldright = [0,0];
    
    for idx = 1:n-1
        [leftpoints(idx, :), rightpoints(idx, :)] = get_border(track.x_X_m(idx), track.y_m(idx), track.x_X_m(idx+1), track.y_m(idx+1), track.w_tr_left_m(idx), track.w_tr_left_m(idx+1), track.w_tr_right_m(idx), track.w_tr_right_m(idx+1), oldleft, oldright);
        oldleft(:) = leftpoints(idx, :);
        oldright(:) = rightpoints(idx, :);
    end
    
    left_border = leftpoints;
    right_border = rightpoints;

end