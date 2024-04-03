close all
clear
clc
carriers = 512;%���ز���
symbols = 14;%OFDM������
c=3*1e8;%����
fc=2.4*1e10;%��Ƶ
delta_f=15000;%���ز����
r=ones(carriers,symbols);
%% ���������ٶ�
range = 117;%m
velocity = 13;%m/s
delta_R=c/(2*carriers*delta_f);%����ֱ���
delta_V=c*delta_f/(2*fc*symbols);%�ٶȷֱ���
%% ����ʱ����Ͷ�������
delay = 2*range/c;
doppler = 2*velocity*fc/c;
k_r=ones(carriers,1);
k_v=ones(1,symbols);
for k=1:carriers
    k_r(k,1)= exp(-1i*2*pi*(k-1)*delta_f*delay);
end
for k=1:symbols
    k_v(1,k)= exp(1i*2*pi*(k-1)*doppler/delta_f); %����ȡ���ز����Ϊ1/T(ofdm)
end
r=r.*(k_r*k_v);
%% ==========Ƶ�׿ն�=====================
matrix_kongdong=zeros(carriers,symbols);
for i=1:symbols
    matrix_kongdong(1:128,i)=r(1:128,i);
    matrix_kongdong(385:512,i)=r(385:512,i);
end
A=ones(carriers,1); 
for i=129:384
    A(i,1)=0;
end
Occupied=256;%ռ�����ز���Ŀ
%% ===========RMS SIMULATION �ٶ�==========
RMSE_R0=zeros(1,31);%Upgrade
RMSE_R1=zeros(1,31);%Traditional
RMSE_R2=zeros(1,31);%CS
K=1000;%�ۼӴ���
for SNR=-30:1:0
    RMSE_0=0;
    RMSE_1=0;
    RMSE_3=0;
    tic
    for k=1:K%�ۼӴ���100��Ҳ���Ը�Ϊ2���ǵú�����Ӧλ���޸�       
       rmatrix_awgn_kongdong=awgn(matrix_kongdong,SNR);%add AWGN
       r_awgn=awgn(r,SNR);
    %% ========Upgrade======
        range_vector = zeros(carriers*10,1);
        r_vector1= zeros(carriers*10,1);
        search_matrix = zeros(carriers*10,1);
        for i=1:symbols
            range_vector(1:carriers,1) = rmatrix_awgn_kongdong(:,i).*A;
            cur_div_IFFT = ifft(range_vector,carriers*10,1);
            search_matrix = search_matrix+abs(cur_div_IFFT);
        end
        search_matrix=search_matrix/symbols;
        [max_R0,index_R0] = max(search_matrix);
        R0=c*(index_R0-1)/(2*carriers*10*delta_f);
     %% ======Traditonal=====
        r_vector1(1:carriers,1)=rmatrix_awgn_kongdong(1:carriers,1);
        cur_div_IFFT1 = ifft(r_vector1,carriers*10,1);
        [max_R1,index_R1] = max(cur_div_IFFT1);
        R1=c*(index_R1-1)/(2*carriers*10*delta_f);

        %����RMSE
        RMSE_0=RMSE_0+(range-R0)*(range-R0);
        RMSE_1=RMSE_1+(range-R1)*(range-R1);

    end
    toc
    RMSE_R0(1,SNR+31)=sqrt(RMSE_0/K);%�����100Ϊǰ���ᵽ��k���ۼӴ����������ۼӴ�������Ӧ���޸ģ�ĿǰΪ2
    RMSE_R1(1,SNR+31)=sqrt(RMSE_1/K); 
end
    SNR_dB=-30:1:0;%������
    figure
    semilogy(SNR_dB,RMSE_R0,'r');
    hold on
    semilogy(SNR_dB,RMSE_R1,'b');
    hold on
    
    
% step1 ����Ĳ�������F���渵��Ҷ�����������ȡ��Ӧ�о�����Ϊ1λ�õ�������      
F=inv(ifft(eye(carriers,carriers))/sqrt(carriers));%������ɢ�渵��Ҷ����
for j=129:384
   F(j,:)=0;
end
F=F(any(F,2),:);%�Ҹ��о����ж�Ӧ���й��ɲ�������A  

%���򻯲�����ȡֵ-30-10dB
lambda=[1,1,1,1,1,1,1,1,9801,8001,8501,7500,6500,6401,5401,5201,6001,6101,4601,5801,...
    3601,5801,5001,6401,5801,5601,5401,5601,5601,5801,5401,5601,5001,4601,5401,5201,5001,5601,5601,5601,5201];
% step2 ������֪��D��div����ÿһ��������Ϊy����Ϊѹ����֪�е�Y
%����ѭ������ÿһ�е�ԭʼ�ź������
    K=50;
    i=1;
    for SNR=-30:1:0
        RMSE_2=0;
        tic
        for k=1:K%�ۼӴ���100��Ҳ���Ը�Ϊ2���ǵú�����Ӧλ���޸�       
          rmatrix_awgn_kongdong=awgn(matrix_kongdong,SNR);%add AWGN
       %% =====����ѹ����֪BPֱ�����ȫƵ��OFDM���ɵľ���-������profile R      Y=FR
            %���
            range_vectorCS = zeros(carriers,1);
            cumulative_R_vectorCS = zeros(carriers,1);
            D1=rmatrix_awgn_kongdong;%D1Ϊ�ع�ǰ��D��div��
            for iii=1:symbols    
                y=D1(:,iii).*A;
                y=y(any(y,2),:);
                x=chat_FISTA(F,y,5000);
                cumulative_R_vectorCS = cumulative_R_vectorCS+abs(x);%���л��۵õ��������
            end
            
            cumulative_R_vectorCS=cumulative_R_vectorCS/symbols;%���ۼӵ��������й�һ��
            [max_R2,index_R2] = max(cumulative_R_vectorCS);
             R2=c*(index_R2-1)/(2*carriers*delta_f);
            RMSE_2=RMSE_2+(range-R2)*(range-R2);
        end
        i=i+1;
        toc
        RMSE_R2(1,SNR+31)=sqrt(RMSE_2/K);        
    end
    SNR_dB=-30:1:0;%������
    semilogy(SNR_dB,RMSE_R2,'g');
    title('Range Measurement Performance');
    xlabel('SNR/dB');
    ylabel('RMSE(r)/m/s');
    legend('Upgrade','Traditonal','FISTA'); 