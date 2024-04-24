function [S_Lr,S_Lt] = getRFChain(AA,L_r,L_t1,G_t,G_r)
        S_Lr = zeros(1,L_r);%选择的接收端射频链路的索引
        S_Lt = zeros(1,L_t1);%选择的发送端射频链路的索引
        w = 1;
        flag = 0;
        while w <= L_r
             location1 = find(abs(AA) == max(max(abs(AA))));
             location2 = location1(1);
             S_Lt(1,w) = ceil(location2/(G_t));
             S_Lr(1,w) = location2-(S_Lt(1,w)-1)*G_t;
              AA(rem((-10 + S_Lr(1,w): 10 + S_Lr(1,w))-1+G_r,G_r)+1,rem((-10+S_Lt(1,w):10+ S_Lt(1,w))-1+G_t,G_t)+1) = 0;
%              AA(rem((--10 + S_Lr(1,w): 5 + S_Lr(1,w))-1+G_r,G_r)+1,rem((-5+S_Lt(1,w):5+ S_Lt(1,w))-1+G_t,G_t)+1) = 0;
             if (w > 1)
                 for w1 = 1:w-1
                      if(S_Lt(1,w) == S_Lt(1,w1))
                          flag = 1;
                          break;
                      elseif(S_Lr(1,w) == S_Lr(1,w1))
                          flag = 1;
                          break;
                      else
                          flag = 0;
                      end
                 end
             end
            w = w + 1;
            if(flag == 1)
                w = w - 1;
            end
        end
end