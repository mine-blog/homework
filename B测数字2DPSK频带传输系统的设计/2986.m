%%%%%%总设计中的频谱波形还没有画出


%%%一些常量
close all;
i=10;  %表明的是10个长度
j=240000;	%对应的是总共的比特数24万比特（10s）
t=linspace(0,10,j);	%表明的是0-10这个区间表现在图像上(时间范围)
fc=4;
fm=i/5;
B=2*fm;

%%%产生基带信号
a=round(rand(1,i));
    for n=1:i 	%对每秒内的信息赋二进制01数
        if a(n)==0	%为0
            for m=j/i*(n-1)+1:j/i*n
                st1(m)=0;	%为0
            end
        else
            for m=j/i*(n-1)+1:j/i*n
                st1(m)=1;	%为1
            end
        end
    end
figure(1);
subplot(311);
plot(t,st1);
title('生成的第一个图像--绝对码');
axis([0,5,-1,2]);
ylabel('绝对码');
xlabel('时间');

%%%差分变换
b=zeros(1,i);	%创造一个10长度的数组
b(1)=a(1);
for n=2:10	%利用的是键控调制中的数字表达式：a(n)=b(n)+b(n-1)
	if a(n)==1
		if b(n-1)==1
			b(n)=0;
		else
			b(n)=1;
		end
	else
		b(n)=b(n-1);
	end
end
st1=t;	%这里将我们指定的空间分为和t相同多的等长矩阵
for n=1:i 	%将得到的值赋给st1矩阵
	if b(n)==0
		for m=j/i*(n-1)+1:j/i*n
			st1(m)=0;
		end
	else
		for m=j/i*(n-1)+1:j/i*n
			st1(m)=1;
		end
	end
end
subplot(312);
plot(t,st1);
title('生成的第二个图像--相对码');
axis([0,5,-1,2]);
ylabel('相对码');
xlabel('时间');

%%%相对码的反码
st2=t;	%再创建一个矩阵
for n=1:j
	if st1(n)==1
		st2(n)=0;
	else
		st2(n)=1;
	end
end
subplot(313)
plot(t,st2)
title('生成的第三个图像--相对码的反码');
axis([0,5,-1,2]);
ylabel('反码');
xlabel('时间');

%%%载波信号
s1=sin(2*pi*fc*t);
figure(2);
plot (s1) ;
title('生成的第四个图像--载波信号');
axis([0,96000,-1,1]);
ylabel('载波信号');
xlabel('时间');

%%%调制信号
d1=st1.*s1;
d2=st2.*(-s1);	%反相180度
figure(3);
subplot(411);
plot(t,d1);
title('生成第五个图像--调制1/2');
ylabel('st1*s1');
xlabel('时间');
axis([0,5,-2,2]);
subplot(412);
plot(t,d2);
title('生成的第六个图像--调制2/2');
ylabel('st2*(-s1)');
xlabel('时间');
axis([0,5,-2,2]);
dpsk=d1+d2;
subplot(413);
plot(t,dpsk);
title('生成的第七个图像--调制后的信号');
ylabel('调制结果');
xlabel('时间');
axis([0,5,-2,2]);
dpsk=awgn(dpsk,10);%加入高斯噪声信噪比是10
subplot(414);
plot(t,dpsk);
title('生成的第八个图像--高斯信道的信号');
xlabel('时间');
ylabel('加噪声的信号');
axis([0,5,-2,2]);
%%%与载波相乘
dpsk=dpsk.*sin(2*pi*fc*t);
figure(4);
subplot(411);
plot(t,dpsk);
title('得到的第九个图像--和载波相乘');
ylabel('和载波相乘');
xlabel('时间');
axis([0,5,-2,2]);

%%%低通滤波
[f,af]=T2F(t,dpsk);	%经过傅里叶变换的信号求频域
%plot(f,af);
title('得到的第十个图像--频域波形');
axis([-20,20,0,3]);
xlabel('频率');
ylabel('幅度');
[t,dpsk]=lpf(f,af,B);	%经过低通滤波过滤掉高频信号
subplot(413);
plot (t,dpsk) ;
title('得到的第十一个图像--低通滤波');
xlabel('时间');
ylabel('通过的低频信号');
axis([0,5,-1,2]);

%%%抽样判决
st=zeros(1,i);
for m=0:i-1
if dpsk(1,m*24000+12000)<0
	st(m+1)=0;
for j=m*24000+1:(m+1)*24000
		dpsk(1,j)=0;
end
else 
	for j=m*24000+1:(m+1)*24000
		st(m+1)=1;
		dpsk(1,j)=1;
	end
end
end
subplot(414);
plot(t,dpsk);
axis([0,5,-1,2]);
title('得到的第十二个图像--抽样判决后波形')
xlabel('时间');
ylabel('判决结果');

%%%码反变换
dt=zeros(1,i);
dt(1)=st(1);
for n=2:i
	if (st(n)-st(n-1))<=0&&(st(n)-st(n-1))>-1
		dt(n)=0;
	else
		dt(n)=1;
	end
end
st=t;
for n=1:i
	if dt(n)<1
		for m=j/i*(n-1)+1:j/i*n
			st(m)=0;
		end
	else
		for m=j/i*(n-1)+1:j/i*n
			st(m)=1;
		end
	end
end
figure(5);
subplot(311);
plot(t,st);
title('解调后的波形');
axis([0,5,-1,2]);

 
function [f,sf]=T2F (t, st) %傅里叶变换
%dt=t(2)-t(1) ;
T=t(end) ;
df=1/T;
N=length(st);
f=-N/2*df:df:N/2*df-df;
sf=fft(st);
sf=T/N*fftshift(sf);
subplot(412);
plot(f,abs(sf(1:round(N))));
%axis([-10,10,0,4]);
title('信号的频谱');
end
 
function [t,st]=F2T(f,sf)   %傅里叶反变换
df=f(2)-f(1);
Fmx=(f(end)-f(1)+df);
dt=1/Fmx;
N=length(sf);
T=dt*N;
t=0:dt:T-dt;
sff=fftshift(sf);
st=Fmx*ifft(sff);
end

function [t,st]=lpf(f,sf,B) %低通滤波器  
df= f (2)-f(1) ;
%T = 1/df ;
hf=zeros(1,length(f)) ;
bf = (-floor( B/df ): floor( B/df )) + floor(length(f)/2 );
hf (bf)=1;
yf=hf.*sf;
[t,st]=F2T(f,yf);
st=real(st);
end


function y=awgn(varargin)
%AWGN Add white Gaussian noise to a signal.
%   Y = AWGN(X,SNR) adds white Gaussian noise to X.  The SNR is in dB.
%   The power of X is assumed to be 0 dBW.  If X is complex, then 
%   AWGN adds complex noise.
%
%   Y = AWGN(X,SNR,SIGPOWER) when SIGPOWER is numeric, it represents 
%   the signal power in dBW. When SIGPOWER is 'measured', AWGN measures
%   the signal power before adding noise.
%
%   Y = AWGN(X,SNR,SIGPOWER,STATE) resets the state of RANDN to STATE.
%
%   Y = AWGN(..., POWERTYPE) specifies the units of SNR and SIGPOWER.
%   POWERTYPE can be 'db' or 'linear'.  If POWERTYPE is 'db', then SNR
%   is measured in dB and SIGPOWER is measured in dBW.  If POWERTYPE is
%   'linear', then SNR is measured as a ratio and SIGPOWER is measured
%   in Watts.
%
%   Example 1: 
%        % To specify the power of X to be 0 dBW and add noise to produce
%        % an SNR of 10dB, use:
%        X = sqrt(2)*sin(0:pi/8:6*pi);
%        Y = awgn(X,10,0);
%
%   Example 2: 
%        % To specify the power of X to be 3 Watts and add noise to
%        % produce a linear SNR of 4, use:
%        X = sqrt(2)*sin(0:pi/8:6*pi);
%        Y = awgn(X,4,3,'linear');
%
%   Example 3: 
%        % To cause AWGN to measure the power of X and add noise to
%        % produce a linear SNR of 4, use:
%        X = sqrt(2)*sin(0:pi/8:6*pi);
%        Y = awgn(X,4,'measured','linear');
%
%   See also WGN, RANDN, and BSC.

%   Copyright 1996-2008 The MathWorks, Inc.
%   $Revision: 1.9.4.6 $  $Date: 2008/08/22 20:23:43 $ 

% --- Initial checks
error(nargchk(2,5,nargin,'struct'));

% --- Value set indicators (used for the string flags)
pModeSet    = 0;
measModeSet = 0;

% --- Set default values
sigPower = 0;
pMode    = 'db';
measMode = 'specify';
state    = [];

% --- Placeholder for the signature string
sigStr = '';

% --- Identify string and numeric arguments
for n=1:nargin
   if(n>1)
      sigStr(size(sigStr,2)+1) = '/';
   end
   % --- Assign the string and numeric flags
   if(ischar(varargin{n}))
      sigStr(size(sigStr,2)+1) = 's';
   elseif(isnumeric(varargin{n}))
      sigStr(size(sigStr,2)+1) = 'n';
   else
      error('comm:awgn:InvalidArg','Only string and numeric arguments are allowed.');
   end
end

% --- Identify parameter signatures and assign values to variables
switch sigStr
   % --- awgn(x, snr)
   case 'n/n'
      sig      = varargin{1};
      reqSNR   = varargin{2};

   % --- awgn(x, snr, sigPower)
   case 'n/n/n'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      sigPower = varargin{3};

   % --- awgn(x, snr, 'measured')
   case 'n/n/s'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      measMode = lower(varargin{3});

      measModeSet = 1;

   % --- awgn(x, snr, sigPower, state)
   case 'n/n/n/n'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      sigPower = varargin{3};
      state    = varargin{4};

   % --- awgn(x, snr, 'measured', state)
   case 'n/n/s/n'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      measMode = lower(varargin{3});
      state    = varargin{4};

      measModeSet = 1;

   % --- awgn(x, snr, sigPower, 'db|linear')
   case 'n/n/n/s'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      sigPower = varargin{3};
      pMode    = lower(varargin{4});

      pModeSet = 1;

   % --- awgn(x, snr, 'measured', 'db|linear')
   case 'n/n/s/s'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      measMode = lower(varargin{3});
      pMode    = lower(varargin{4});

      measModeSet = 1;
      pModeSet    = 1;

   % --- awgn(x, snr, sigPower, state, 'db|linear')
   case 'n/n/n/n/s'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      sigPower = varargin{3};
      state    = varargin{4};
      pMode    = lower(varargin{5});

      pModeSet = 1;

   % --- awgn(x, snr, 'measured', state, 'db|linear')
   case 'n/n/s/n/s'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      measMode = lower(varargin{3});
      state    = varargin{4};
      pMode    = lower(varargin{5});

      measModeSet = 1;
      pModeSet    = 1;

   otherwise
      error('comm:awgn:InvalidSyntax','Syntax error.');
end   

% --- Parameters have all been set, either to their defaults or by the values passed in,
%     so perform range and type checks

% --- sig
if(isempty(sig))
   error('comm:awgn:NoInput','An input signal must be given.');
end

if(ndims(sig)>2)
   error('comm:awgn:InvalidSignalDims','The input signal must have 2 or fewer dimensions.');
end

% --- measMode
if(measModeSet)
   if(~strcmp(measMode,'measured'))
      error('comm:awgn:InvalidSigPower','The signal power parameter must be numeric or ''measured''.');
   end
end

% --- pMode
if(pModeSet)
   switch pMode
   case {'db' 'linear'}
   otherwise
      error('comm:awgn:InvalidPowerType','The signal power mode must be ''db'' or ''linear''.');
   end
end

% -- reqSNR
if(any([~isreal(reqSNR) (length(reqSNR)>1) (isempty(reqSNR))]))
   error('comm:awgn:InvalidSNR','The signal-to-noise ratio must be a real scalar.');
end

if(strcmp(pMode,'linear'))
   if(reqSNR<=0)
      error('comm:awgn:InvalidSNRForLinearMode','In linear mode, the signal-to-noise ratio must be > 0.');
   end
end

% --- sigPower
if(~strcmp(measMode,'measured'))

   % --- If measMode is not 'measured', then the signal power must be specified
   if(any([~isreal(sigPower) (length(sigPower)>1) (isempty(sigPower))]))
      error('comm:awgn:InvalidSigPower','The signal power value must be a real scalar.');
   end
   
   if(strcmp(pMode,'linear'))
      if(sigPower<0)
         error('comm:awgn:InvalidSigPowerForLinearMode','In linear mode, the signal power must be >= 0.');
      end
   end

end

% --- state
if(~isempty(state))
   if(any([~isreal(state) (length(state)>1) (isempty(state)) any((state-floor(state))~=0)]))
      error('comm:awgn:InvaildState','The State must be a real, integer scalar.');
   end
end

% --- All parameters are valid, so no extra checking is required

% --- Check the signal power.  This needs to consider power measurements on matrices
if(strcmp(measMode,'measured'))
   sigPower = sum(abs(sig(:)).^2)/length(sig(:));

   if(strcmp(pMode,'db'))
      sigPower = 10*log10(sigPower);
   end
end

% --- Compute the required noise power
switch lower(pMode)
   case 'linear'
      noisePower = sigPower/reqSNR;
   case 'db'
      noisePower = sigPower-reqSNR;
      pMode = 'dbw';
end

% --- Add the noise
if(isreal(sig))
   opType = 'real';
else
   opType = 'complex';
end

y = sig+wgn(size(sig,1), size(sig,2), noisePower, 1, state, pMode, opType);

end