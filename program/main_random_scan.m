clear variables;
close all;

% parameters
r = 0.3; % 既知の画素値の割合
p = 5; % p x p の画像パッチをモデル化する。奇数でないといけない。

A = zeros(5,35);

l = 35;
m = 5;

for s = 1:l

    for t = 1:m

        % 画像の読み込みとグレースケールへの変換、サイズ取得
        I_color_true = imread('cloud.png');
        I_gray_true = double(rgb2gray(I_color_true));
        [I_x, I_y] = size(I_gray_true);
        
        % 既知の画素の位置を表す行列Omの生成（ランダムに生成）
        Om = make_Om_random_scan(I_x,I_y,r);
        
        % 欠損画像の生成
        I_ms = I_gray_true.*Om;
        
        % 画像の表示
        %figure
        %imshow(uint8(I_gray_true))
        %title('元の画像')
        %figure
        %imshow(uint8(I_ms))
        %title('欠損画像')
        
        
        % 既知の画素を表す行列Omから、ベクトルbの既知の要素を表現するベクトルObを生成
        bias = (p+1)/2;
        Ob = reshape( Om(bias:I_x-bias+1,bias:I_y-bias+1)',[],1 );
        
        % ベクトルbの未知の要素を表現するベクトルOcを生成
        Oc = ~Ob; % ベクトルObの０と１を反転させたもの
        
        
        % データ行列Xとベクトルbを画像から生成し、初期値として保存
        [X0, b0] = make_X_b(I_ms,p);
        
        
        %-------------------------------
        % ARモデルに基づく画像修復アルゴリズム
        % 現在は、周囲の画素の平均値を代入するだけなので、これを改良する
        %-------------------------------
        I_result = I_ms; % 修復画像の初期値は欠損画像
        
        %a = ones(p*p-1,1)/(p*p-1); % 周りの画素値の平均値を計算する係数ベクトル

        n = s;
        
        for u = 1:n
            % データ行列Xとベクトルbを画像I_resultから生成
            [X , b] = make_X_b(I_result,p);
        
            % Xとbからaを計算
            a = X\b;
        
            % 係数ベクトルaを使ってbを計算
            b = X*a;
        
            % ベクトルbの未知の要素の値はそのままで、既知の要素の値を代入
            b = b.*Oc + b0.*Ob;
        
            % ベクトルbから画像I_msを生成
            I_result(bias:I_x-bias+1,bias:I_y-bias+1) = reshape(b, I_y-(bias-1)*2,I_x-(bias-1)*2)';
        
            
        end
    
        % 画素あたりの誤差計算（ただし、周辺部を除く）
        error = norm( I_result(bias:I_x-bias+1,bias:I_y-bias+1)-I_gray_true(bias:I_x-bias+1,bias:I_y-bias+1),'fro')/(I_x*I_y);
    
        A(t,s) = error;

    end

end

B = mean(A);

%-------------------------------


% 結果表示
%figure
%imshow(uint8(I_result))
%title('修復結果')



