function [bestModel, maxGoodCount] = ransac_fundamental_matrix(m1,m2,threshold,niters,confidence)

%% 利用Ransac 方法计算基本矩阵

maxGoodCount=0;
    for testNo=1:niters
        ind = randperm(size(m1,1),8);

        F = estimate_fundamental_matrix(m1(ind,:), m2(ind,:));
        err = computererror(m1, m2, F);
        goodCount = 0; %%初始化内点个数为零
        for i=1:size(m1,1)
            if err(i) < threshold^2
                goodCount = goodCount+1;
            end
        end

        if goodCount> max(maxGoodCount,8-1)
            bestModel = F;
            maxGoodCount = goodCount;
            %% 更新迭代次数
            niters = RANSACUpdateNumIters(confidence, (size(m1,1) - goodCount)/size(m1,1), 8);
        end

    end
end

%% 该函数用来衰减迭代次数
function niters = RANSACUpdateNumIters(p,ep,modelPoints)
niters = round(log(1-p)/log(1-ep^modelPoints));
end


function err = computererror(m1, m2, F)
    for i=1:size(m1,1)
        left = m1(i,:);
        right = m2(i,:);
        left(:,3)=1;
        right(:,3)=1;

        a = left*F(1,:)';
        b = left*F(2,:)';
        c = left*F(3,:)';
        d2 = [a,b,c]*(right');
        s2 = 1/(a*a+b*b);

        a = right*F(:,1);
        b = right*F(:,2);
        c = right*F(:,3);
        d1 = [a,b,c]*(left');
        s1 = 1/(a*a+b*b);
        err(i) = max(d1*d1*s1, d2*d2*s2);
    end
end
