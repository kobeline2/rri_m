%%

T1 = 10;
T2 = 5;
[T3,T1] = Func1(T1, T2)


function [T3,T1] = Func1(T1,T2)
    T4 = T1 + T2;
    T1 = T1 - 2;
    T3 = T4 + 10;
end
