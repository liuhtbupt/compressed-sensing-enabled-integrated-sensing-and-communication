close all
clear
clc
tic
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
rmatrix_awgn_kongdong=add_awgn_noise(matrix_kongdong,10);

    %% =========traditional ==========
     r_vt=zeros(carriers,1);
     search_vt=zeros(carriers,1);
     for i=1:symbols
         r_vt=rmatrix_awgn_kongdong(:,i);
         r_ifft=ifft(r_vt);
         search_vt=search_vt+abs(r_ifft);
     end
     search_vt=search_vt/symbols;
     [max_R1,index_R1] = max(search_vt);
     R1=c*(index_R1-1)/(2*carriers*delta_f);
    %% =========Upgrade=============
    range_vector = zeros(carriers,1);
      search_matrix = zeros(carriers,1);
        for i=1:symbols
            range_vector = rmatrix_awgn_kongdong(:,i).*A;
            cur_div_IFFT = ifft(range_vector);
            search_matrix = search_matrix+abs(cur_div_IFFT);
        end
        search_matrix=search_matrix/symbols;
        [max_R0,index_R0] = max(search_matrix);
        R0=c*(index_R0-1)/(2*carriers*delta_f);
  %% =====FISTA==========================================
      addpath('utils/');
         D1=rmatrix_awgn_kongdong;%D1Ϊ�ع�ǰ��D��div��
      % ���ò�����Ԫ�ص�ѹ������B������ѹ���ع����õ��ع�����D2
        %step1 ������ɢFourier����ȡ�õ���������A���������ز�������256����ȫ�����ز�������512����
           F=ifft(eye(carriers,carriers))/sqrt(carriers);%������ɢ�渵��Ҷ����
           F=inv(F);%����ɢ�渵��Ҷ����������
           for j=129:384
               F(j,:)=0;
           end
           F;%�Ҹ��о����ж�Ӧ���й��ɲ�������A           
        % step2 ��y=ceita*x��֪����������ѹ������о���y���ع���ԭʼ�о���x
            %����ѭ������ÿһ�е�ԭʼ�ź������
            D2=zeros(carriers,symbols);
            for iii=1:symbols
                y=D1(:,iii).*A;             
                 x=chat_FISTA(F,y,5201);
               D2(:,iii)=x;
            end
          
              %���
               range_vectorCS = zeros(carriers,1);
               cumulative_R_vectorCS = zeros(carriers,1);
               for i=1:symbols
                    range_vectorCS = D2(:,i);
                    cur_div_IFFTCS = range_vectorCS;
                    cumulative_R_vectorCS = cumulative_R_vectorCS+abs(cur_div_IFFTCS);
               end
               cumulative_R_vectorCS=cumulative_R_vectorCS/symbols;%���ۼӵ��������й�һ��
               [max_R3,index_R3] = max(cumulative_R_vectorCS);
               R3=c*(index_R3-1)/(2*carriers*delta_f);
              
               %% =====ѹ����֪====
        figure(1)
        y1=abs(cumulative_R_vectorCS/max_R3);
        y2=abs(search_matrix/max_R0);
        y3=abs(search_vt/max_R1);
        plot(y1,'r');
        hold on
        plot(y2,'b');
        hold on
        plot(y3,'y');
        xlabel('index_R');
        ylabel('Energy spectrum');
        title('NC-OFDM range measurement');
        legend('����FISTA�㷨','�Ľ���2D FFT�㷨','��ͳ��2D FFT')
  toc     
%     zp=BaseZoom();
%              zp.plot;
%  %% ======= ��ȫOFDM���й����׹��Ƶõ��ľ���=======
%   x=zeros(512,14);
%   sdsd=zeros(512,1);
%   for i= 1:14
%       sd=r(:,i);
%       sd_ifft=ifft(sd);
%       x(:,i)=sd_ifft;
%   end
%   [x_,y_] = meshgrid(1:14,1:512);
%   mesh(x_,y_,abs(x))
%   xlswrite('��ʵx����.xlsx',real(x),1,'A2:N513')
%   xlswrite('��ʵx����.xlsx',imag(x),1,'O2:AB513')
%   
%   D10=matrix_kongdong;
%   for i=129:384
%       D10(i,:)=0;
%   end 
%   aaa=chat_FISTA(F,D10(:,1),400);
%   plot(abs(aaa))
%   %%
%   [best,lam,sc] = range_Kfold(F,D10,x);
%   best
%   figure(1)
%   plot(lam,sc)
%   plot(abs(x(:,1)))
%   %%
%   D10=matrix_kongdong;
%   for i=129:384
%       D10(i,:)=0;
%   end
%   xlswrite('range_mode1_û�������ն�����.xlsx',real(D10),1,'A2:N513')
%   xlswrite('range_mode1_û�������ն�����.xlsx',imag(D10),1,'O2:AB513')
%     xlswrite('range_mode1_��ɢ����Ҷ����.xlsx',real(F),1,'A2:SR513')
%   xlswrite('range_mode1_��ɢ����Ҷ����.xlsx',imag(F),1,'SS2:AMJ513')