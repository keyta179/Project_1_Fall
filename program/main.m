clear variables;
close all;

% parameters
r = 0.3; % ���m�̉�f�l�̊���
p = 5; % p x p �̉摜�p�b�`�����f��������B��łȂ��Ƃ����Ȃ��B

% �摜�̓ǂݍ��݂ƃO���[�X�P�[���ւ̕ϊ��A�T�C�Y�擾
I_color_true = imread('cloud.png');
% I_color_true_R = I_color_true(:,:,1);
% I_color_true_G = I_color_true(:,:,2);
% I_color_true_B = I_color_true(:,:,3);
I_gray_true = double(rgb2gray(I_color_true));
[I_x, I_y] = size(I_gray_true);

% ���m�̉�f�̈ʒu��\���s��Om�̐����i�����_���ɐ����j
Om = make_Om(I_x,I_y,r);

% �����摜�̐���
I_ms = I_gray_true.*Om;

% �����摜�̐����A�e�F�̏��擾(�J���[)
% I_ms_c_R = I_color_true_R.*uint8(Om);
% I_ms_c_G = I_color_true_G.*uint8(Om);
% I_ms_c_B = I_color_true_B.*uint8(Om);
% I_ms_c(:,:,1) = I_ms_c_R;
% I_ms_c(:,:,2) = I_ms_c_G;
% I_ms_c(:,:,3) = I_ms_c_B;

% �摜�̕\��
figure
imshow(uint8(I_gray_true))
title('���̉摜')
figure
imshow(uint8(I_ms))
title('�����摜')

% �摜�̕\��(�J���[)
% figure
% imshow(uint8(I_color_true))
% title('���̉摜')
% figure
% imshow(uint8(I_ms_c))
% title('�����摜')


% ���m�̉�f��\���s��Om����A�x�N�g��b�̊��m�̗v�f��\������x�N�g��Ob�𐶐�
bias = (p+1)/2;
Ob = reshape( Om(bias:I_x-bias+1,bias:I_y-bias+1)',[],1 );

% �x�N�g��b�̖��m�̗v�f��\������x�N�g��Oc�𐶐�
Oc = ~Ob; % �x�N�g��Ob�̂O�ƂP�𔽓]����������


% �f�[�^�s��X�ƃx�N�g��b���摜���琶�����A�����l�Ƃ��ĕۑ�
[X0, b0] = make_X_b(I_ms,p);

% �f�[�^�s��X�ƃx�N�g��b���摜���琶�����A�����l�Ƃ��ĕۑ�(�J���[)
% [X0_R, b0_R] = make_X_b(I_ms_c_R,p);
% [X0_G, b0_G] = make_X_b(I_ms_c_G,p);
% [X0_B, b0_B] = make_X_b(I_ms_c_B,p);


%-------------------------------
% AR���f���Ɋ�Â��摜�C���A���S���Y��
% ���݂́A���͂̉�f�̕��ϒl�������邾���Ȃ̂ŁA��������ǂ���
%-------------------------------
I_result = I_ms; % �C���摜�̏����l�͌����摜
% a = ones(p*p-1,1)/(p*p-1); % ����̉�f�l�̕��ϒl���v�Z����W���x�N�g��

% I_result_R = I_ms_c_R;
% I_result_G = I_ms_c_G;
% I_result_B = I_ms_c_B;

% �������v�Z����s��𐶐�

% D = make_D(I_ms,p);

% Dt = transpose(D);

% [x_I, y_I] = size(D);

% I = eye(x_I,x_I);

% lambda = 1;

% D_result = I + lambda*D*Dt;

n = 50;

for k = 1:n
    % �f�[�^�s��X�ƃx�N�g��b���摜I_result���琶��
    [X , b] = make_X_b(I_result,p);
   
    % �f�[�^�s��X�ƃx�N�g��b���摜I_result���琶��(�J���[)

    % [X_R , b_R] = make_X_b(I_result_R,p);
    % [X_G , b_G] = make_X_b(I_result_G,p);
    % [X_B , b_B] = make_X_b(I_result_B,p);

    % X��b����a���v�Z
    a = X\b;

    % X��b����a���v�Z(�J���[)
    % a_R = X_R\b_R;
    % a_G = X_G\b_G;
    % a_B = X_B\b_B;

    % �W���x�N�g��a���g����b���v�Z
    b = X*a;

    % b = D_result\X*a;

    % �W���x�N�g��a���g����b���v�Z(�J���[)
    % b_R = X_R*a_R;
    % b_G = X_G*a_G;
    % b_B = X_B*a_B;

    % �x�N�g��b�̖��m�̗v�f�̒l�͂��̂܂܂ŁA���m�̗v�f�̒l����
    b = b.*Oc + b0.*Ob;

    % �x�N�g��b�̖��m�̗v�f�̒l�͂��̂܂܂ŁA���m�̗v�f�̒l����(�J���[)
    % b_R = b_R.*Oc + b0_R.*Ob;
    % b_G = b_G.*Oc + b0_G.*Ob;
    % b_B = b_B.*Oc + b0_B.*Ob;

    % �x�N�g��b����摜I_ms�𐶐�
    I_result(bias:I_x-bias+1,bias:I_y-bias+1) = reshape(b, I_y-(bias-1)*2,I_x-(bias-1)*2)';

    % �x�N�g��b����摜I_ms�𐶐�(�J���[)
    % I_result_R(bias:I_x-bias+1,bias:I_y-bias+1) = reshape(b_R, I_y-(bias-1)*2,I_x-(bias-1)*2)';
    % I_result_G(bias:I_x-bias+1,bias:I_y-bias+1) = reshape(b_G, I_y-(bias-1)*2,I_x-(bias-1)*2)';
    % I_result_B(bias:I_x-bias+1,bias:I_y-bias+1) = reshape(b_B, I_y-(bias-1)*2,I_x-(bias-1)*2)';
end

% �ĂуJ���[�摜���擾
% I_result_c(:,:,1) = I_result_R;
% I_result_c(:,:,2) = I_result_G;
% I_result_c(:,:,3) = I_result_B;

%-------------------------------


% ���ʕ\��
figure
imshow(uint8(I_result))
title('�C������')

% ��f������̌덷�v�Z�i�������A���ӕ��������j
error = norm( I_result(bias:I_x-bias+1,bias:I_y-bias+1)-I_gray_true(bias:I_x-bias+1,bias:I_y-bias+1),'fro')/(I_x*I_y);

% ���ʕ\��(�J���[)
% figure
% imshow(uint8(I_result_c))
% title('�C������')
