%% 多くの引数を関数にわたす際の方法
% a, b, c, d, eという引数を渡す関数を想定
% 関数は全て最下段で定義（MATLABのルール）
% 全て1e7回関数を回して、計算時間を比較する。

%% 方法1
% 普通に引数を渡す
a = 1;
b = 2;
c = 3;
d = 4;
e = 5;
f = 6;

tic
for I = 1:1e7
    out2 = func2(C);
end
toc

%% 方法2
% 構造体にして渡す（見通しがよい）

C.a = 1;
C.b = 2;
C.c = 3;
C.d = 4;
C.e = 5;
C.f = 6;

tic
for I = 1:1e7
    out3 = func3();
end
toc

%% 方法3
% グローバル変数にしてしまう

global ag bg cg dg eg fg
ag = 1;
bg = 2;
cg = 3;
dg = 4;
eg = 5;
fg = 6;

tic 
for I = 1:1e7
    out1 = func1(a, b, c, d, e, f);
end
toc

%%
% 方法1
function out = func1(a, b, c, d, e, f)
    out = a * sin(b) * cos(c) / d^e + f;
end

% 方法2
function out = func2(C)
    out = C.a * sin(C.b) * cos(C.c) / C.d^C.e + C.f;
end

% 方法3
function out = func3()
    global ag bg cg dg eg fg
    out = ag * sin(bg) * cos(cg) / dg^eg + fg;
end