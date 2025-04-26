# application-week-6-p0-solved
**TO GET THIS SOLUTION VISIT:** [Application Week 6 P0 Solved](https://www.ankitcodinghub.com/product/application-p0-solved-6/)


---

📩 **If you need this solution or have special requests:** **Email:** ankitcoding@gmail.com  
📱 **WhatsApp:** +1 419 877 7882  
📄 **Get a quote instantly using this form:** [Ask Homework Questions](https://www.ankitcodinghub.com/services/ask-homework-questions/)

*We deliver fast, professional, and affordable academic help.*

---

<h2>Description</h2>



<div class="kk-star-ratings kksr-auto kksr-align-center kksr-valign-top" data-payload="{&quot;align&quot;:&quot;center&quot;,&quot;id&quot;:&quot;124730&quot;,&quot;slug&quot;:&quot;default&quot;,&quot;valign&quot;:&quot;top&quot;,&quot;ignore&quot;:&quot;&quot;,&quot;reference&quot;:&quot;auto&quot;,&quot;class&quot;:&quot;&quot;,&quot;count&quot;:&quot;3&quot;,&quot;legendonly&quot;:&quot;&quot;,&quot;readonly&quot;:&quot;&quot;,&quot;score&quot;:&quot;5&quot;,&quot;starsonly&quot;:&quot;&quot;,&quot;best&quot;:&quot;5&quot;,&quot;gap&quot;:&quot;4&quot;,&quot;greet&quot;:&quot;Rate this product&quot;,&quot;legend&quot;:&quot;5\/5 - (3 votes)&quot;,&quot;size&quot;:&quot;24&quot;,&quot;title&quot;:&quot;Application Week 6 P0 Solved&quot;,&quot;width&quot;:&quot;138&quot;,&quot;_legend&quot;:&quot;{score}\/{best} - ({count} {votes})&quot;,&quot;font_factor&quot;:&quot;1.25&quot;}">

<div class="kksr-stars">

<div class="kksr-stars-inactive">
            <div class="kksr-star" data-star="1" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" data-star="2" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" data-star="3" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" data-star="4" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" data-star="5" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
    </div>

<div class="kksr-stars-active" style="width: 138px;">
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
    </div>
</div>


<div class="kksr-legend" style="font-size: 19.2px;">
            5/5 - (3 votes)    </div>
    </div>
Application of Matlab for Finance

Class 6

Today’s Class

I Simulation

I Normal Stock Price Modelhttps://.com I

Log Normal Stock Price Model

I Black-Scholes Option Pricing

Add

Random Number Generators

I rand(m,n): uniform random number on the interval (0,1)

I randi(m,n): uniform integer random number on the interval (0,1) I randn(m,n): standard Normal random number

I normrnd(mu,sigma,m,n)https://.com: normal random number with mean mu and standard deviation sigma

I trnd(nu,m,n): student t-distribution random number with nu

degrees of freedomAdd I randg(m,n): standard Gamma random number

I (m,n) defined the output matrix size, m-by-n, that stores the simulated numbers.

Examples: Normal &amp; Student-t

I Simulate 100 × 1 standard Normal random variables

I Simulate 100 × 1 Normal random variables from N(0.08,0.22) use both randn and normrnd

X N(µ,σ2)

https://.comX = µ + σ N(0,1)

I Simulate 100 × 1 random variables from a Student t with 8 degrees of freedom

1 % Simulate 100 number from N(0,1) x = randn(100,1); 23 % Simulate 100 number from N(0.08, 4) % X

4 = mu + sigma* N(0,1) y = 0.08 + 0.2 *

5 randn(100,1); y2 =

6 normrnd(0.08,0.2,[100,1]);

7 % Simulate student t−distribution x = trnd(8,100,1);

8

9

10

Add

Set Seeds for Random Generator

I Sometimes, we want to use the same sequence of random numbers for various reasons, such as code debugging or to generate

reproductive results.

I We can set the seed for the random number generator usinghttps://.com rng.

I Call the stored seed every time when we want to regenerate the same sequence of random numbers.

1 % seeding the random number generator s = rng; % set

2 the seed for generator a = rand(1,5) rng(s); % call the

3 stored seed b = rand(1,5)

4

5

6

Add

Simulate Asset Prices

I In finance, the price of a particular stock at a future time t is unknown

at the present

I We often think of it as being a random variable St

I If we simplly assumehttps://.comSt follows a normal distribution N(µ,σ) for t = 1,2,3,.., then

S Add S

I In expectation, E0 since both 1 and 2 are random draws from N(0,1).

I We need the time dimension variations.

Normal Stock Price Model

I We assume the stock price follows a stochastic process

I Thenhttps://.comSt

t

I This is process known as the Brownian Motion

I Then gross return and net return on stock as

Add poSt w√ coder

Rt = = 1 + µ∆t + σ ∆t

S0

Rtnet

Simulate the Stock Price Process

I Consider stock with annual return of 0.15 and annual volatility of Assignment 0.3, today you observe its priceP2501 roject Exam Help$1. Simulate a 1 year path for this stock price with ∆t = year.

1

2

4 dt = 1/250;

5

6 T = 1; tgrid = 0: dt: T;

7 N = length(tgrid);

8

9 % Set up parameters &amp; initialize the price vector S0 = 1; mu = 0.15;

10

11

12

13

14

15

17 sigma = 0.3;

18

19 S = zeros(1,N);

% Simulate random number epsilon eps = randn(1,N);

% Simulate stock prices

S = S0*(1 + mu*tgrid + sigma*sqrt(tgrid).*eps); plot(tgrid, S) legend( ‘S’) xlabel(‘Time(yr)’) ylabel(‘Asset Price($)’)

Simulate the Stock Price Process: 3 Stocks

13 S0 = 1;

14

15 % Simulate 3 random numbers epsilon

16 eps = randn(3,N);

18 % Simulate stock prices

19 S1 = S0 *(1+ mu1*tgrid + sigma1*sqrt(tgrid).* eps(1,:));

20 S2 = S0 *(1+ mu2*tgrid + sigma2*sqrt(tgrid).* eps(2,:));

21 S3 = S0 *(1+ mu3*tgrid + sigma3*sqrt(tgrid).* eps(3,:));

22 plot(tgrid, S1, ‘−’, tgrid, S2,’:’, tgrid, S3,’−−’)

23 legend( ‘S1’, ‘S2’, ‘S3’)

24 xlabel(‘Time(yr)’)

25 ylabel(‘Asset Price($)’)

Normal Stock Price Model

https://.com Add

Log-Normal Stock Price Model

I Issues with Normal stock price: negative stock prices

I wherein

, and ∆t = t − 0 = t

I St follows a log-normal distribution, as often referred as Geometric

Brownian

Motion:https://.com Add

St St

rt = ln(St) − ln(S0) = ln( )

S0

Log Normal Stock Price Model

4

5 % dt = 1/250 yr, T = 1 yr,

6 7

9 dt = 1/250;

10

11 % time gride = T/dt + 1 = 251 grid points

12 T = 1; tgrid = 0: dt: T;

13 N = length(tgrid);

14 S1 = zeros(1,N); S2 = zeros(1,N); S3 = zeros(1,N);= 1; S0

15

16 % Simulate random number epsilon eps =

17 randn(3,N);

18

19

20

21

S1=S0Add

*exp((mu1−0.5*(sigma1^2))*tgrid+sigma1*sqrt(tgrid).*eps); S2=S0*exp((mu2−0.5*(sigma2^2))*tgrid+sigma2*sqrt(tgrid).*eps); S3=S0*exp((mu3−0.5*(sigma3^2))*tgrid+sigma3*sqrt(tgrid).*eps);

plot(tgrid, S1, ‘−’, tgrid,

S2,’:’, tgrid, S3,’−−’) legend( ‘S1’, ‘S2’, ‘S3’) xlabel(‘Time(yr)’) ylabel(‘Asset Price($)’)

Log-Normal Stock Price Model

https://.com Add

Options Pricing

I V(S,t) is the value of an option

I C(S,t): Call options give the right to purchase the underlying asset at

I S the value of stock price (i.e. underlying asset)

I Khttps://.comthe strike price of the option contract (i.e. the agreed price)

I rAdd is the risk free rate

I σ is the volatility of the underlying stock.

C(S,T) = max(S −K,0) P(S,T) = max(K −S,0)

Option Pricing Exercises Option Pricing Simulation

I Use the normal stock price model simulate 10,000,000 scenarios for a stock

that with S0 = 100,T = 1,µ = 10%,σ = 20%;

I Calculate the expected price of an European Call and Put option on the this stock with K = 100,r = 5%

https://.comI Today’s price is the discounted value of expected future

payoffs

I Law of Large Number implies the sample mean converts to the true mean of the underlying distribution.

I Note: in the previous exercises, we simulate over a sequence of time grids in the future.Add

1

Option Pricing Exercises

1 S0 = 100; % Value of the underlying

2 K = 100; % Strike (exercise price)

3 T = 1; % Maturity

4 mu = 0.10; % Stock price mean

5

6

7 sigma =

10

13

14

15

16

r = 0.10; % Risk free interest rate

% Number of simulation trials

eps = randn(M,1); %

*exp((mu−0.5*(sigma^2))*T+sigma*sqrt(T).*eps);

S_T = S0*(1 + mu*T+ sigma*sqrt(T).*eps);

payoff_call=max(S_T−K,0); −S_T,0);

Option Pricing Exercises

p_call = mean(exp(−r*T)*payoff_call); %

continuous discount p_put =

mean(exp(−r*T)*payoff_put);

Options Price 2: The Black-Scholes Formula

I The price of a call option is given by

C(T) = SN(d1) −Ke−rTN(d2)

I The price of a put option is given by https://.comP(T) = Ke−rTN(−d2) −SN(−d1)

I where ln(S/K) + (r + σ2/2)T d1 = √ Add σ pT owcoder

√

Option Pricing Exercises

d2 = d1 − σ T

I and N(d1) and N(d2) denotes the standard cumulative normal probability for the values of d1 and d2. It is the probability that a random draw from a normal distribution.

I Exercise 2:

I Create a function perform the Black-Scholes Formula to determine the above options’ price

https://.comI The function shall allows the user to define a option as

call or put

I Compare with the above method

I Currently, the normal stock price process is used in the simulation based method;

Add I The difference among the two comes from the different underlying I Black-Scholes model assumes a log-normal stock price; stock price processes;

Option Pricing Exercises

2: Function

1 function price = BlackScholesPrice(S, K, T, r, sigma, CallorPut) 2 % This function calculates option price base on the Black−Schole formula. 3

4 % Input: S: spot stock price

5 % K: strike price

6 % T: maturity

7 % r: interest rate

string input as ‘Call’ or ‘Put’ option

% sigma: volatility

10

11 if strcmp(CallorPut,’Call’) == 1

12 phi = 1;

14 phi =−1;

Option Pricing Exercises

15 else

error(‘Invalid Option Type’)

18

19 d1 = (log(S/K) + (r + 0.5 * sigma^2)* T)./ sigma.* sqrt(T);

20 Nd1 = normcdf(phi*d1,0,1); 21

22 d2 = d1− sigma.* sqrt(T); 23

Nd2 = normcdf(phi*d2,0,1); 24

25 price = phi.*S.*Nd1− phi.*K.*exp(−r * T).*Nd2;

26 end

2: Main Command

1 %%%%%%%%%% Option parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

2 S = 100; % Value of the underlying

3

4

ption Pricing Exercises5

6 K = 100;

7

8 (T = 1;exercise price ) % Maturity r = 0.10; % Risk free 9 % Strike

10 interest rate sigma = 0.20; % Volatility %

11 parameters = [S, K, T, r, sigma];

12

% call functions to calculate the price of the option

CallPrice = BlackScholesPrice(S, K, T, r, sigma, ‘Call’);

PutPrice = BlackScholesPrice(S, K, T, r, sigma, ‘Put’);

% fprintf output fprintf(‘The Price of a European Call is : %.2f ‘, CallPrice);

Add

Take Away

I Basic Simulation Code

I Simulate stock prices follow normal and log normal process

I Pricing Options using simulation and Black-Scholes Modelhttps://.com

I In your coursework, you will use the log normal stock price for simulation. Will there be any difference between the option price based on the simulation versus the Black-Scholes?Add
