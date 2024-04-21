function [Om] = make_Om_random_line_scan(m,n,r)


% ランダムラインスキャン



k = 400; % ランダムラインスキャンの回数

Om = zeros(m,n);
for s = 1:k
    a = rand(1);
    c = randsample(2,1);
    if c == 1
        a = -a;
    end
    b = randperm(m,1);
    d = randsample(2,1);
    if d == 1
        b = -b;
    end
    for y = 1:n
        x = int16(a*y + b);
        if x > 0 && x <= m
            Om(x,y) = 1;
        end
    end
end