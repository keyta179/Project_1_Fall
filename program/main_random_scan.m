clear variables;
close all;

% parameters
r = 0.3; % ���m�̉�f�l�̊���
p = 5; % p x p �̉摜�p�b�`�����f��������B��łȂ��Ƃ����Ȃ��B

A = zeros(5,35);

l = 35;
m = 5;

for s = 1:l

    for t = 1:m

        % �摜�̓ǂݍ��݂ƃO���[�X�P�[���ւ̕ϊ��A�T�C�Y�擾
        I_color_true = imread('cloud.png');
        I_gray_true = double(rgb2gray(I_color_true));
        [I_x, I_y] = size(I_gray_true);
        
        % ���m�̉�f�̈ʒu��\���s��Om�̐����i�����_���ɐ����j
        Om = make_Om_random_scan(I_x,I_y,r);
        
        % �����摜�̐���
        I_ms = I_gray_true.*Om;
        
        % �摜�̕\��
        %figure
        %imshow(uint8(I_gray_true))
        %title('���̉摜')
        %figure
        %imshow(uint8(I_ms))
        %title('�����摜')
        
        
        % ���m�̉�f��\���s��Om����A�x�N�g��b�̊��m�̗v�f��\������x�N�g��Ob�𐶐�
        bias = (p+1)/2;
        Ob = reshape( Om(bias:I_x-bias+1,bias:I_y-bias+1)',[],1 );
        
        % �x�N�g��b�̖��m�̗v�f��\������x�N�g��Oc�𐶐�
        Oc = ~Ob; % �x�N�g��Ob�̂O�ƂP�𔽓]����������
        
        
        % �f�[�^�s��X�ƃx�N�g��b���摜���琶�����A�����l�Ƃ��ĕۑ�
        [X0, b0] = make_X_b(I_ms,p);
        
        
        %-------------------------------
        % AR���f���Ɋ�Â��摜�C���A���S���Y��
        % ���݂́A���͂̉�f�̕��ϒl�������邾���Ȃ̂ŁA��������ǂ���
        %-------------------------------
        I_result = I_ms; % �C���摜�̏����l�͌����摜
        
        %a = ones(p*p-1,1)/(p*p-1); % ����̉�f�l�̕��ϒl���v�Z����W���x�N�g��

        n = s;
        
        for u = 1:n
            % �f�[�^�s��X�ƃx�N�g��b���摜I_result���琶��
            [X , b] = make_X_b(I_result,p);
        
            % X��b����a���v�Z
            a = X\b;
        
            % �W���x�N�g��a���g����b���v�Z
            b = X*a;
        
            % �x�N�g��b�̖��m�̗v�f�̒l�͂��̂܂܂ŁA���m�̗v�f�̒l����
            b = b.*Oc + b0.*Ob;
        
            % �x�N�g��b����摜I_ms�𐶐�
            I_result(bias:I_x-bias+1,bias:I_y-bias+1) = reshape(b, I_y-(bias-1)*2,I_x-(bias-1)*2)';
        
            
        end
    
        % ��f������̌덷�v�Z�i�������A���ӕ��������j
        error = norm( I_result(bias:I_x-bias+1,bias:I_y-bias+1)-I_gray_true(bias:I_x-bias+1,bias:I_y-bias+1),'fro')/(I_x*I_y);
    
        A(t,s) = error;

    end

end

B = mean(A);

%-------------------------------


% ���ʕ\��
%figure
%imshow(uint8(I_result))
%title('�C������')



