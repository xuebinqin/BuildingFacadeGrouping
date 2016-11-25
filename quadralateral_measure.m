function qm = quadralateral_measure(p1,p2,p3,p4)
% input for points and return the quadralateral measurement
% p_i are all 1-by-2 vector   
    pred_p4 = (p2-p1)+(p3-p1)+p1;
%     qm = norm((p4-pred_p4),2)/(norm((p2-p1)+(p3-p1),2));
    qm = norm((p4-pred_p4),2);
end