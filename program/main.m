clear variables;
close all;

% parameters
r = 0.3; % 既知の画素値の割合
p = 5; % p x p の画像パッチをモデル化する。奇数でないといけない。

% 画像の読み込みとグレースケールへの変換、サイズ取得
I_color_true = imread('cloud.png');
% I_color_true_R = I_color_true(:,:,1);
% I_color_true_G = I_color_true(:,:,2);
% I_color_true_B = I_color_true(:,:,3);
I_gray_true = double(rgb2gray(I_color_true));
[I_x, I_y] = size(I_gray_true);

% 既知の画素の位置を表す行列Omの生成（ランダムに生成）
Om = make_Om(I_x,I_y,r);

% 欠損画像の生成
I_ms = I_gray_true.*Om;

% 欠損画像の生成、各色の情報取得(カラー)
% I_ms_c_R = I_color_true_R.*uint8(Om);
% I_ms_c_G = I_color_true_G.*uint8(Om);
% I_ms_c_B = I_color_true_B.*uint8(Om);
% I_ms_c(:,:,1) = I_ms_c_R;
% I_ms_c(:,:,2) = I_ms_c_G;
% I_ms_c(:,:,3) = I_ms_c_B;

% 画像の表示
figure
imshow(uint8(I_gray_true))
title('元の画像')
figure
imshow(uint8(I_ms))
title('欠損画像')

% 画像の表示(カラー)
% figure
% imshow(uint8(I_color_true))
% title('元の画像')
% figure
% imshow(uint8(I_ms_c))
% title('欠損画像')


% 既知の画素を表す行列Omから、ベクトルbの既知の要素を表現するベクトルObを生成
bias = (p+1)/2;
Ob = reshape( Om(bias:I_x-bias+1,bias:I_y-bias+1)',[],1 );

% ベクトルbの未知の要素を表現するベクトルOcを生成
Oc = ~Ob; % ベクトルObの０と１を反転させたもの


% データ行列Xとベクトルbを画像から生成し、初期値として保存
[X0, b0] = make_X_b(I_ms,p);

% データ行列Xとベクトルbを画像から生成し、初期値として保存(カラー)
% [X0_R, b0_R] = make_X_b(I_ms_c_R,p);
% [X0_G, b0_G] = make_X_b(I_ms_c_G,p);
% [X0_B, b0_B] = make_X_b(I_ms_c_B,p);


%-------------------------------
% ARモデルに基づく画像修復アルゴリズム
% 現在は、周囲の画素の平均値を代入するだけなので、これを改良する
%-------------------------------
I_result = I_ms; % 修復画像の初期値は欠損画像
% a = ones(p*p-1,1)/(p*p-1); % 周りの画素値の平均値を計算する係数ベクトル

% I_result_R = I_ms_c_R;
% I_result_G = I_ms_c_G;
% I_result_B = I_ms_c_B;

% 差分を計算する行列を生成

% D = make_D(I_ms,p);

% Dt = transpose(D);

% [x_I, y_I] = size(D);

% I = eye(x_I,x_I);

% lambda = 1;

% D_result = I + lambda*D*Dt;

n = 50;

for k = 1:n
    % データ行列Xとベクトルbを画像I_resultから生成
    [X , b] = make_X_b(I_result,p);
   
    % データ行列Xとベクトルbを画像I_resultから生成(カラー)

    % [X_R , b_R] = make_X_b(I_result_R,p);
    % [X_G , b_G] = make_X_b(I_result_G,p);
    % [X_B , b_B] = make_X_b(I_result_B,p);

    % Xとbからaを計算
    a = X\b;

    % Xとbからaを計算(カラー)
    % a_R = X_R\b_R;
    % a_G = X_G\b_G;
    % a_B = X_B\b_B;

    % 係数ベクトルaを使ってbを計算
    b = X*a;

    % b = D_result\X*a;

    % 係数ベクトルaを使ってbを計算(カラー)
    % b_R = X_R*a_R;
    % b_G = X_G*a_G;
    % b_B = X_B*a_B;

    % ベクトルbの未知の要素の値はそのままで、既知の要素の値を代入
    b = b.*Oc + b0.*Ob;

    % ベクトルbの未知の要素の値はそのままで、既知の要素の値を代入(カラー)
    % b_R = b_R.*Oc + b0_R.*Ob;
    % b_G = b_G.*Oc + b0_G.*Ob;
    % b_B = b_B.*Oc + b0_B.*Ob;

    % ベクトルbから画像I_msを生成
    I_result(bias:I_x-bias+1,bias:I_y-bias+1) = reshape(b, I_y-(bias-1)*2,I_x-(bias-1)*2)';

    % ベクトルbから画像I_msを生成(カラー)
    % I_result_R(bias:I_x-bias+1,bias:I_y-bias+1) = reshape(b_R, I_y-(bias-1)*2,I_x-(bias-1)*2)';
    % I_result_G(bias:I_x-bias+1,bias:I_y-bias+1) = reshape(b_G, I_y-(bias-1)*2,I_x-(bias-1)*2)';
    % I_result_B(bias:I_x-bias+1,bias:I_y-bias+1) = reshape(b_B, I_y-(bias-1)*2,I_x-(bias-1)*2)';
end

% 再びカラー画像を取得
% I_result_c(:,:,1) = I_result_R;
% I_result_c(:,:,2) = I_result_G;
% I_result_c(:,:,3) = I_result_B;

%-------------------------------


% 結果表示
figure
imshow(uint8(I_result))
title('修復結果')

% 画素あたりの誤差計算（ただし、周辺部を除く）
error = norm( I_result(bias:I_x-bias+1,bias:I_y-bias+1)-I_gray_true(bias:I_x-bias+1,bias:I_y-bias+1),'fro')/(I_x*I_y);

% 結果表示(カラー)
% figure
% imshow(uint8(I_result_c))
% title('修復結果')
