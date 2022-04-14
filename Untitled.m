testIntKMLoc;
% 构建矩阵,包含所有可能的约束点
v = 1:12;
pu = nchoosek(v,6);
% 循环,找出来不独立的约束点集
for i = 1: length(pu)
if rank(K(pu(i,:),pu(i,:)))== length(K(pu(i,:),pu(i,:)))
jug(i,:) = 1;
else
jug(i,:) = 0;
end
end
errorindex = pu(jug==0,:);
% 去除只约束同方向的约束点集
count = 1;
for i = 1:length(errorindex)
    if (length(setxor(errorindex(i,:),[1 4 7 10]))==6 && ...
        length(setxor(errorindex(i,:),[2 5 8 11]))==6 && ...
        length(setxor(errorindex(i,:),[3 6 9 12]))==6)
        temp(count,:) = errorindex(i,:);
        count = count + 1;
    end
end

%% 检测约束点集pu1是否独立
pu1=[1 3 5];
disp('秩为');
disp(rank(K(pu1,pu1)))
disp('方阵大小为');
disp(size(K(pu1,pu1),1))