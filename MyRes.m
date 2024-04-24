function [Res,s] = MyRes(pinv_RR_W,y_w,RR_W)
    s = pinv_RR_W*y_w;
    Res = y_w-RR_W*s;
 end