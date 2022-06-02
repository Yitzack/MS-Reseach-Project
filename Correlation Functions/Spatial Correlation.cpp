#include<iostream>
#include<fstream>
#include<cmath>
using namespace std;

pair<long double,long double> Integrate(long double, long double(*)(long double, long double));
long double Kernel(long double, long double, long double);

long double Interpolation(long double, long double, long double[151][259]);	//(s,P) comes in, f(s,P) comes out
long double i_index(long double, long double);				//Turns (s,P) into index variables
long double j_index(long double, long double);
long double Basis0(long double);						//BSpline Basis functions
long double Basis1(long double);
long double Basis2(long double);
long double Basis3(long double);
long double Basisn(long double);

long double Interacting(long double, long double);		//Returns the Interacting Spectral Function
long double Non_interacting(long double, long double);	//Returns the Non-interacting Spectral Function

void Import(char*);

//Global varibles. I would rather not do this. But it needs to be done to ensure that Interacting() and Non_interacting() have the same function pointer shape
long double ImGf[151][259], ImGvCf[151][259], ImGvLf[151][259], ImGvQf[151][259], ImGVf[151][259];	//Imaginary contributions to the spectral function
long double ReGvCf[151][259], ReGvLf[151][259], ReGvQf[151][259], ReGVf[151][259];			//Real contributions to the spectral function
long double Coupling_Fraction, Vacuum_Coupling;							//Fraction of Vacuum Coupling Counstant, Vacuum Coupling Constant

int main(int argc, char* argv[])
{
	pair<long double,long double> holder;
	Import(argv[1]);

	for(long double z = .25; z <= 6.25; z+=.25)
	{
		holder = Integrate(z, Non_interacting);
		cout << "{" << z << ",Around[" << holder.first << "," << holder.second << "],Around[" << flush;
		holder = Integrate(z, Interacting);
		cout << holder.first << "," << holder.second << "]}" << endl;
	}

	return(0);
}

pair<long double,long double> Integrate(long double z, long double(*sigma)(long double, long double))
{
#if ORDER == 7
	long double Disp[] = {0.5773502691896257645, 0.9258200997725514616};	//Displacement from center for 16th order Gauss-Kronrod integration from Mathematica
	long double w[] = {28./45.,27./55.,98./495.};	//Weight of the function at Disp
	long double werr[] = {28./45.,-28./55.,98./495.};	//Weight of the error estimate function at Disp
#elif ORDER == 16
	long double Disp[] = {0.27963041316178319341346652274898, 0.53846931010568309103631442070021, 0.75416672657084922044081716694612, 0.90617984593866399279762687829939, 0.98408536009484246449617293463614};	//Displacement from center for 16th order Gauss-Kronrod integration from Mathematica
	long double w[] = {0.28298741785749121320425560137111, 0.272849801912558922340993264484456, 0.24104033922864758669994261122326, 0.186800796556492657467800026878486, 0.115233316622473394024626845880574, 0.0425820367510818328645094508476701};	//Weight of the function at Disp
	long double werr[] = {0.272849801912558922340993264484456, -0.28590147103139767568463328751778, 0.272849801912558922340993264484456, -0.237588331270718881341348903612376, 0.186800796556492657467800026878486, -0.121693568433715693489637194839344, 0.0425820367510818328645094508476701};	//Weight of the error estimate function at Disp
#elif ORDER == 37
	long double Disp[] = {0.12523340851146891547244136946385, 0.24850574832046927626779096036272, 0.36783149899818019375269153664372, 0.48133945047815709293594361501883, 0.58731795428661744729670241894053, 0.68405989547005589394492910034115, 0.76990267419430468703689383321282, 0.84355812416115324479214188505984, 0.90411725637047485667846586611910, 0.95053779594312129654906019513162, 0.98156063424671925069054909014928, 0.99693392252959542691235023725839};	//Displacement from center for 37th order Gauss-Kronrod integration
	long double w[] = {0.125556893905474335304296132860078, 0.12458416453615607343731247320923, 0.121626303523948383246099758091310, 0.11671205350175682629358074530573, 0.110022604977644072635907398742250, 0.101649732279060277715688770491228, 0.0915494682950492105281719397396142, 0.079920275333601701493392609529783, 0.0672509070508399303049409400473161, 0.053697017607756251228889163320458, 0.0389152304692994771150896322858629, 0.023036084038982232591084580367969, 0.00825771143316839575769392243921158};	//Weights for 37th order Gauss-Kronrod integration
	long double werr[] = {0.125556893905474335304296132860078, -0.12456288127724671156324996283372, 0.121626303523948383246099758091310, -0.116780483036597982467269153619148, 0.110022604977644072635907398742250, -0.101517694444005644033375685318571, 0.0915494682950492105281719397396142, -0.080158053209744524841259920013576, 0.0672509070508399303049409400473161, -0.053242308387562179731365554873538, 0.0389152304692994771150896322858629, -0.0241392523475295946035313811170481, 0.00825771143316839575769392243921158};	//Error weights for 37th order Gauss-Kronrod integration
#elif ORDER == 64
	long double Disp[] = {0.07296724759579058908806031209955, 0.14556185416089509093703098233869, 0.21738009387172883298179561608565, 0.28802131680240109660079251606460, 0.35712525365992107730319197431219, 0.42434212020743878357366888854379, 0.48929590889566864140105140228151, 0.55161883588721980705901879672431, 0.61099712378976603905186076209981, 0.66713880419741231930596666999034, 0.71972345257202475786195418640521, 0.76843996347567790861587785130623, 0.81305258084753191324907332122919, 0.85336336458331728364725063858757, 0.88912831533118433426270707847375, 0.92009933415040082879018713371497, 0.94614721677851262706535437885119, 0.96722683856630629431662221490770, 0.98318486356168595383392464989066, 0.99375217062038950026024203593794, 0.99896193068731718686198704308270};	//Displacement from center for 64th order Gauss-Kronrod integration
	long double w[] = {0.07302748796276908297314068625699, 0.0728456528305084410294911786777899, 0.072275482222034062651328264004603, 0.0712933463637845561021040045941832, 0.069929765008398753379866909680647, 0.0682207062610252150704081432085748, 0.066149577815951110013778554529288, 0.0636956152491007998422517557993921, 0.060898680204661962686034547272808, 0.0578104574976095223026959977828257, 0.054418675178383649262446811973774, 0.0506971977965175119748827957014245, 0.046697324731514497924831858406713, 0.0424969328811184855040644730480614, 0.038083153453465002583875378860838, 0.0334025463208076220711022803409528, 0.028518426570204953082884284481435, 0.0235741990830290305646428964274167, 0.018560055408876216262191994606715, 0.0133046426070883174480910406380580, 0.0078183373021241792478434900552220, 0.00279548123241156950861299678078312};	//Weights for 64th order Gauss-Kronrod integration
	long double werr[] = {-0.073053645686921344218844461426381, 0.0728456528305084410294911786777899, -0.072248921767935996412498902549150, 0.0712933463637845561021040045941832, -0.069957629782674401342266514186936, 0.0682207062610252150704081432085748, -0.066119360817386351767274019967487, 0.0636956152491007998422517557993921, -0.060932735849066571509332629852926, 0.0578104574976095223026959977828257, -0.054378623988764728401027766096332, 0.0506971977965175119748827957014245, -0.0467470987245193636284578827072191, 0.0424969328811184855040644730480614, -0.0380169601749142994331762744393456, 0.0334025463208076220711022803409528, -0.0286159988566522552007515419910130, 0.0235741990830290305646428964274167, -0.0183937343619762775377586736926149, 0.0133046426070883174480910406380580, -0.0081988909556501540763811268032490, 0.00279548123241156950861299678078312};	//Error weights for 64th order Gauss-Kronrod integration*/
#elif ORDER == 97
	long double Disp[] = {0.04830766568773831623481257044050, 0.09650269687689436580083125301636, 0.14447196158279649348518637359881, 0.19210360898314249727164159225259, 0.23928736225213707454460320916550, 0.28591245858945975941660710119035, 0.33186860228212764977991680573019, 0.37704942115412110544533548554981, 0.42135127613063534536411943617243, 0.46466930848199221775617820037971, 0.50689990893222939002374747437782, 0.54794631419915247868093950160974, 0.58771575724076232904074547640183, 0.62611293770182399782023837975795, 0.66304426693021520097511516866324, 0.69842655779521049288477014244748, 0.73218211874028968038742666509127, 0.76422825199780370415066012687437, 0.79448379596794240696309729897043, 0.82288295013605132164826884947213, 0.84936761373256997013369300496774, 0.87386976894531060612966180261363, 0.89632115576605212396530724371921, 0.91667726665136432427534565885792, 0.93490607593773968917091913483541, 0.95095468484866118538988275061533, 0.96476225558750643077381192811827, 0.97631028361466380719766964319312, 0.98561151154526833540017504463090, 0.99262803526297191268579115643847, 0.99726386184948156354498112866504, 0.99954590212436447863561028016112};	//Displacement from center for 97th order Gauss-Kronrod integration
	long double w[] = {0.0483263839865677583754454340005196, 0.04827019307577738559871207883299, 0.0481009691854577469278465438391581, 0.047818908736988472212263584657302, 0.0474260618738823823628799498539221, 0.046922968281703611103480713352515, 0.0463087567380257132403812984724432, 0.045585826564547070280575461963265, 0.0447586387497669372951991920753785, 0.043827544030139749046816152673798, 0.0427911155964467469336549254437386, 0.041654019985643051398296413297779, 0.0404234923703730966723492694215763, 0.039099420133306611207482127044236, 0.0376791306456133985148959737488490, 0.036169769475642299860958391303290, 0.0345821227447330341307263834168846, 0.032915077643903600263296476785556, 0.0311633255619737371711558486849664, 0.029336956689620661368615614585584, 0.0274520984222104037831477067959546, 0.025505695480894652814528899899308, 0.0234866596721633245920879128647720, 0.021408913184821915955777518020231, 0.0192987714303268112944037397635239, 0.017149805209784253256085826220233, 0.0149361036060860273850967513325094, 0.012676054806654402859368883399502, 0.0104239873988068188280342507620158, 0.0081725040385316684143438051192845, 0.00584173707916669330394797663999954, 0.0034268187757723709355745755785249, 0.00122336081795147180029303715064886};	//Weights for 97th order Gauss-Kronrod integration
	long double werr[] = {0.0483263839865677583754454340005196, -0.04826989543895041496805275123058, 0.0481009691854577469278465438391581, -0.047819811342286387206818417546829, 0.0474260618738823823628799498539221, -0.046921430799100954535699524315602, 0.0463087567380257132403812984724432, -0.045588052131216814432293115148372, 0.0447586387497669372951991920753785, -0.043824548974264062095955310078004, 0.0427911155964467469336549254437386, -0.041657904241303703823902661306570, 0.0404234923703730966723492694215763, -0.039094475653763695264258791784071, 0.0376791306456133985148959737488490, -0.036176024633206206364440965175198, 0.0345821227447330341307263834168846, -0.032907145132458246574353586921383, 0.0311633255619737371711558486849664, -0.029347136788914885776668022714587, 0.0274520984222104037831477067959546, -0.025492363781481523381634344790214, 0.0234866596721633245920879128647720, -0.0214269848374047647011011285858944, 0.0192987714303268112944037397635239, -0.0171240577032371798466019060321402, 0.0149361036060860273850967513325094, -0.0127160105026076565963837063897218, 0.0104239873988068188280342507620158, -0.0081018906923740021908267570871021, 0.00584173707916669330394797663999954, -0.0035917912336977256648324881603283, 0.00122336081795147180029303715064886};	//Error weights for 97th order Gauss-Kronrod integration
#endif

	long double e1, e2, P1, P2, c1, c2;
	long double P_interval = M_PI/(4.*z);
	long double a = 0, b = .1, c, d;
	long double P_Max;
	long double Answer = 0, Partial_Answer;
	long double Error = 0, Partial_Error;
	long double holder;
	int i, j, k;

	if(sigma == Interacting)
		P_Max = int(32.*z/M_PI)*M_PI/z;
	else
		P_Max = int(70.4*z/M_PI)*M_PI/z;

	do
	{
		c = 0;
		d = .1;
		do	//General parallelogram
		{
			Partial_Answer = 0;
			Partial_Error = 0;
#if ORDER == 7
			for(i = 0; i < 2; i++) //Integrate the sub-interval
#elif ORDER == 16
			for(i = 0; i < 5; i++)
#elif ORDER == 37
			for(i = 0; i < 12; i++)
#elif ORDER == 64
			for(i = 0; i < 21; i++)
#elif ORDER == 97
			for(i = 0; i < 32; i++)
#endif
			{
				c1 = (b+a-Disp[i]*(b-a))/2.;
				c2 = (b+a+Disp[i]*(b-a))/2.;
#if ORDER == 7
				for(j = 0; j < 2; j++) //Integrate the sub-interval
#elif ORDER == 16
				for(j = 0; j < 5; j++)
#elif ORDER == 37
				for(j = 0; j < 12; j++)
#elif ORDER == 64
				for(j = 0; j < 21; j++)
#elif ORDER == 97
				for(j = 0; j < 32; j++)
#endif
				{
					e1 = (c+d-Disp[j]*(d-c))/2.;
					e2 = (c+d+Disp[j]*(d-c))/2.;

					P1 = e1+c1;	//First for c1
					P2 = e2+c1;
					holder = sigma(e1, P1)*Kernel(e1, P1, z)*(b-a)*(d-c)/4.;
					Partial_Answer += w[i]*w[j]*holder;
					Partial_Error += werr[i]*werr[j]*holder;
					holder = sigma(e2, P2)*Kernel(e2, P2, z)*(b-a)*(d-c)/4.;
					Partial_Answer += w[i]*w[j]*holder;
					Partial_Error += werr[i]*werr[j]*holder;

					P1 = e1+c2;	//Second for c1
					P2 = e2+c2;
					holder = sigma(e1, P1)*Kernel(e1, P1, z)*(b-a)*(d-c)/4.;
					Partial_Answer += w[i]*w[j]*holder;
					Partial_Error += werr[i]*werr[j]*holder;
					holder = sigma(e2, P2)*Kernel(e2, P2, z)*(b-a)*(d-c)/4.;
					Partial_Answer += w[i]*w[j]*holder;
					Partial_Error += werr[i]*werr[j]*holder;
				}

				P1 = (c+d)/2.+c1;	//First for c1
				P2 = (c+d)/2.+c2;	//Second for c2
				holder = sigma((c+d)/2., P1)*Kernel((c+d)/2., P1, z)*(b-a)*(d-c)/4.;	//First for c1
				Partial_Answer += w[i]*w[j]*holder;
				Partial_Error += werr[i]*werr[j]*holder;

				holder = sigma((c+d)/2., P2)*Kernel((c+d)/2., P2, z)*(b-a)*(d-c)/4.;	//Second for c2
				Partial_Answer += w[i]*w[j]*holder;
				Partial_Error += werr[i]*werr[j]*holder;

			}

#if ORDER == 7
			for(j = 0; j < 2; j++) //Integrate the sub-interval
#elif ORDER == 16
			for(j = 0; j < 5; j++)
#elif ORDER == 37
			for(j = 0; j < 12; j++)
#elif ORDER == 64
			for(j = 0; j < 21; j++)
#elif ORDER == 97
			for(j = 0; j < 32; j++)
#endif
			{
				e1 = (c+d-Disp[j]*(d-c))/2.;
				e2 = (c+d+Disp[j]*(d-c))/2.;

				P1 = e1+(a+b)/2.;	//Second for c1
				P2 = e2+(a+b)/2.;
				holder = sigma(e1, P1)*Kernel(e1, P1, z)*(b-a)*(d-c)/4.;
				Partial_Answer += w[i]*w[j]*holder;
				Partial_Error += werr[i]*werr[j]*holder;
				holder = sigma(e2, P2)*Kernel(e2, P2, z)*(b-a)*(d-c)/4.;
				Partial_Answer += w[i]*w[j]*holder;
				Partial_Error += werr[i]*werr[j]*holder;
			}
			P1 = (a+b)/2.+(c+d)/2.;
			holder = sigma((c+d)/2., P1)*Kernel((c+d)/2., P1, z)*(b-a)*(d-c)/4.;
			Partial_Answer += w[i]*w[j]*holder;
			Partial_Error += werr[i]*werr[j]*holder;
			
			c = d;
			if(d < 20.8)
				d += .1;
			else
				d += 1.;

			if(d > P_Max-b)
				d = P_Max-b;

			Answer += Partial_Answer;
			Error += abs(Partial_Error);
		}while(c+b < P_Max);

		Partial_Answer = 0;	//One last triangle to make the trapazoid
		Partial_Error = 0;
		c = P_Max-b;
		d = P_Max-a;
#if ORDER == 7
		for(i = 0; i < 2; i++) //Integrate the sub-interval
#elif ORDER == 16
		for(i = 0; i < 5; i++)
#elif ORDER == 37
		for(i = 0; i < 12; i++)
#elif ORDER == 64
		for(i = 0; i < 21; i++)
#elif ORDER == 97
		for(i = 0; i < 32; i++)
#endif
		{
			c1 = (b+a-Disp[i]*(b-a))/2.;
			c2 = (b+a+Disp[i]*(b-a))/2.;
#if ORDER == 7
			for(j = 0; j < 2; j++) //Integrate the sub-interval
#elif ORDER == 16
			for(j = 0; j < 5; j++)
#elif ORDER == 37
			for(j = 0; j < 12; j++)
#elif ORDER == 64
			for(j = 0; j < 21; j++)
#elif ORDER == 97
			for(j = 0; j < 32; j++)
#endif
			{
				e1 = (c+P_Max-c1-Disp[j]*(P_Max-c1-c))/2.;	//First for c1
				e2 = (c+P_Max-c1+Disp[j]*(P_Max-c1-c))/2.;
				P1 = e1+c1;
				P2 = e2+c1;
				holder = sigma(e1, P1)*Kernel(e1, P1, z)*(b-a)*(P_Max-c1-c)/4.;
				Partial_Answer += w[i]*w[j]*holder;
				Partial_Error += werr[i]*werr[j]*holder;
				holder = sigma(e2, P2)*Kernel(e2, P2, z)*(b-a)*(P_Max-c1-c)/4.;
				Partial_Answer += w[i]*w[j]*holder;
				Partial_Error += werr[i]*werr[j]*holder;

				e1 = (c+P_Max-c2-Disp[j]*(P_Max-c2-c))/2.;	//Second for c1
				e2 = (c+P_Max-c2+Disp[j]*(P_Max-c2-c))/2.;
				P1 = e1+c2;
				P2 = e2+c2;
				holder = sigma(e1, P1)*Kernel(e1, P1, z)*(b-a)*(P_Max-c2-c)/4.;
				Partial_Answer += w[i]*w[j]*holder;
				Partial_Error += werr[i]*werr[j]*holder;
				holder = sigma(e2, P2)*Kernel(e2, P2, z)*(b-a)*(P_Max-c2-c)/4.;
				Partial_Answer += w[i]*w[j]*holder;
				Partial_Error += werr[i]*werr[j]*holder;
			}

			e1 = (c+P_Max-c1)/2.;	//First for c1
			e2 = (c+P_Max-c2)/2.;	//Second for c2
			P1 = e1+c1;
			P2 = e2+c2;
			holder = sigma(e1, P1)*Kernel(e1, P1, z)*(b-a)*(P_Max-c1-c)/4.;	//First for c1
			Partial_Answer += w[i]*w[j]*holder;
			Partial_Error += werr[i]*werr[j]*holder;

			holder = sigma(e2, P2)*Kernel(e2, P2, z)*(b-a)*(P_Max-c2-c)/4.;	//Second for c2
			Partial_Answer += w[i]*w[j]*holder;
			Partial_Error += werr[i]*werr[j]*holder;
		}
#if ORDER == 7
		for(j = 0; j < 2; j++) //Integrate the sub-interval
#elif ORDER == 16
		for(j = 0; j < 5; j++)
#elif ORDER == 37
		for(j = 0; j < 12; j++)
#elif ORDER == 64
		for(j = 0; j < 21; j++)
#elif ORDER == 97
		for(j = 0; j < 32; j++)
#endif
		{
			e1 = (c+P_Max-(a+b)/2.-Disp[j]*(P_Max-(a+b)/2.-c))/2.;
			e2 = (c+P_Max-(a+b)/2.+Disp[j]*(P_Max-(a+b)/2.-c))/2.;

			P1 = e1+(a+b)/2.;	//Second for c1
			P2 = e2+(a+b)/2.;
			holder = sigma(e1, P1)*Kernel(e1, P1, z)*(b-a)*(P_Max-(a+b)/2.-c)/4.;
			Partial_Answer += w[i]*w[j]*holder;
			Partial_Error += werr[i]*werr[j]*holder;
			holder = sigma(e2, P2)*Kernel(e2, P2, z)*(b-a)*(P_Max-(a+b)/2.-c)/4.;
			Partial_Answer += w[i]*w[j]*holder;
			Partial_Error += werr[i]*werr[j]*holder;
		}
		P1 = (a+b)/2.+(c+P_Max-(a+b)/2.)/2.;
		holder = sigma((c+P_Max-(a+b)/2.)/2., P1)*Kernel((c+d)/2., P1, z)*(b-a)*(P_Max-(a+b)/2.-c)/4.;
		Partial_Answer += w[i]*w[j]*holder;
		Partial_Error += werr[i]*werr[j]*holder;

		Answer += Partial_Answer;
		Error += abs(Partial_Error);

		a = b;
		b += .1;
		if(b > 15)
			b = 15;
	}while(b < 15 && a != b);

	return(pair<long double, long double>(Answer, Error));
}

long double Kernel(long double e, long double P, long double z)
{
	return(2.*cos(P*z)/e);
}

long double Interacting(long double e, long double P)	//Returns the Interacting Spectral Function
{
	long double ImGvC = Interpolation(e,P,ImGvCf);	//Get the interpolation of the contributions
	long double ImGvL = sqrt(1.5)*Interpolation(e,P,ImGvCf);
	long double ImGvQ = Interpolation(e,P,ImGvCf)/2.;
	long double ImGV = Interpolation(e,P,ImGVf);
	long double ReGvC = Interpolation(e,P,ReGvCf);
	long double ReGvL = Interpolation(e,P,ReGvLf);
	long double ReGvQ = Interpolation(e,P,ReGvQf);
	long double ReGV = Interpolation(e,P,ReGVf);

	//Each of the 3 interacting spectral function contributions
	long double sigmaC = 3./M_PI/4.*Coupling_Fraction*Vacuum_Coupling*(ImGV*pow(ReGvC,2)-pow(ImGvC,2)*ImGV-2.*(ReGV-1.)*ImGvC*ReGvC)/(pow(ReGV-1.,2)+pow(ImGV,2));
	long double sigmaL = 3./M_PI*Coupling_Fraction*Vacuum_Coupling*(ImGV*pow(ReGvL,2)-pow(ImGvL,2)*ImGV-2.*(ReGV-1.)*ImGvL*ReGvL)/(pow(ReGV-1.,2)+pow(ImGV,2));
	long double sigmaQ = 3./M_PI*Coupling_Fraction*Vacuum_Coupling*(ImGV*pow(ReGvQ,2)-pow(ImGvQ,2)*ImGV-2.*(ReGV-1.)*ImGvQ*ReGvQ)/(pow(ReGV-1.,2)+pow(ImGV,2));

	return(sigmaC+sigmaL+sigmaQ);
}

long double Non_interacting(long double e, long double P)	//Returns the Non-interacting Spectral Function
{
	return(-3./M_PI*Interpolation(e, P, ImGf));
}

long double Interpolation(long double e, long double P, long double f[151][259])
{
	long double i = i_index(e,P);					//Index i
	long double j = j_index(e,P);					//Index j
	if(i < 0 || i > 150 || j < 0 || j > 258)	//Prevent bad data returns
	{
		cerr << e << " " << pow(e,2)-pow(P,2) << " " << P << " " << i_index(e,P) << " " << j_index(e,P) << endl;
		return(0./0.);
	}

	int offset_i = -1, offset_j = -1;					 //Index offsets in calling up control points, normally -1, but can be 0 or -2 on the ends
	long double (*Basisi[4])(long double) = {Basisn,Basisn,Basisn,Basisn}; //Basis Functions in the i direction
	long double (*Basisj[4])(long double) = {Basisn,Basisn,Basisn,Basisn}; //Basis Functions in the j direction
	long double zx[4][4], zy[4][4];					 //z from the x direction and y direction
	long double answer = 0;						 //Result

	if(i < 2)	//Reassign the function pointers for the x/i direction
	{
		Basisi[0] = Basis0;
		Basisi[1] = Basis1;
		Basisi[2] = Basis2;
		Basisi[3] = Basis3;
		if(i < 1)
			offset_i = 0;
	}
	else if(i < 3)
	{
		Basisi[0] = Basis1;
		Basisi[1] = Basis2;
		Basisi[2] = Basis3;
	}
	else if(i < 4)
	{
		Basisi[0] = Basis2;
		Basisi[1] = Basis3;
	}
	else if(i < 5)
	{
		Basisi[0] = Basis3;
	}
	else if(148 <= i)
	{
		Basisi[0] = Basis3;
		Basisi[1] = Basis2;
		Basisi[2] = Basis1;
		Basisi[3] = Basis0;
		if(149 <= i)
			offset_i = -2;
	}
	else if(147 <= i)
	{
		Basisi[1] = Basis3;
		Basisi[2] = Basis2;
		Basisi[3] = Basis1;
	}
	else if(146 <= i)
	{
		Basisi[2] = Basis3;
		Basisi[3] = Basis2;
	}
	else if(145 <= i)
	{
		Basisi[3] = Basis3;
	}

	if(j < 2 || (j >= 208 && j <= 210))	//Reassign the function pointers for the y/j direction
	{
		Basisj[0] = Basis0;
		Basisj[1] = Basis1;
		Basisj[2] = Basis2;
		Basisj[3] = Basis3;
		if(j < 1 || (j >= 208 && j < 209))
			offset_j = 0;
	}
	else if(j < 3 || (j >= 210 && j < 211))
	{
		Basisj[0] = Basis1;
		Basisj[1] = Basis2;
		Basisj[2] = Basis3;
	}
	else if(j < 4 || (j >= 211 && j < 212))
	{
		Basisj[0] = Basis2;
		Basisj[1] = Basis3;
	}
	else if(j < 5 || (j >= 212 && j < 213))
	{
		Basisj[0] = Basis3;
	}
	else if((206 <= j && j < 208) || 256 <= j)
	{
		Basisj[0] = Basis3;
		Basisj[1] = Basis2;
		Basisj[2] = Basis1;
		Basisj[3] = Basis0;
		if((207 <= j && j < 208) || 257 <= j)
			offset_j = -2;
	}
	else if((205 <= j && j < 208) || 255 <= j)
	{
		Basisj[1] = Basis3;
		Basisj[2] = Basis2;
		Basisj[3] = Basis1;
	}
	else if((204 <= j && j < 208) || 254 <= j)
	{
		Basisj[2] = Basis3;
		Basisj[3] = Basis2;
	}
	else if((203 <= j && j < 208) || 253 <= j)
	{
		Basisj[3] = Basis3;
	}

	for(int i_count = 0; i_count < 4; i_count++)	//Evaluate the Basis Functions
		for(int j_count = 0; j_count < 4; j_count++)
		{
			if(Basisj[j_count] != Basisn && j < 5)
				zy[i_count][j_count] = Basisj[j_count](j);
			else if(Basisj[j_count] != Basisn && (208 <= j && j < 213))
				zy[i_count][j_count] = Basisj[j_count](j-208);
			else if(Basisj[j_count] != Basisn && (203 <= j && j < 208))
				zy[i_count][j_count] = Basisj[j_count](208-j);
			else if(Basisj[j_count] != Basisn && 253 <= j)
				zy[i_count][j_count] = Basisj[j_count](258-j);
			else
				zy[i_count][j_count] = Basisj[j_count](i-int(i)+j_count+2.);	//Some how the fractional part of the arguments were trading dimensions. Dispite trying to be careful and keeping x/i in the first argument and y/j in the second argument, I was ending up with f(int(i)+j-int(j),int(j)+i-int(i)) instead of f(int(i)+i-int(i),int(j)+j-int(j))=f(i,j)

			if(Basisi[i_count] != Basisn && i < 5)
				zx[i_count][j_count] = Basisi[i_count](i);
			else if(Basisi[i_count] != Basisn && 145 <= i)
				zx[i_count][j_count] = Basisi[i_count](150-i);
			else
				zx[i_count][j_count] = Basisi[i_count](j-int(j)+i_count+2.);

			answer += zx[i_count][j_count]*zy[i_count][j_count]*f[int(i)+i_count+offset_i][int(j)+j_count+offset_j];
		}

	return(answer);
}

long double j_index(long double e, long double P)
{
	if(e <= 20.8)
		return(10.*e);
	else
		return(187.2+e);
}

long double i_index(long double e, long double P)
{
	return(10.*(P-e));
}

long double Basis0(long double x)
{
	if(0 <= x && x <= 2)
		return(-pow((x-2.)/2.,3));
	return(0);
}

long double Basis1(long double x)
{
	if(0 <= x && x < 2)
		return(x*(19.*pow(x,2)-90.*x+108.)/72.);
	else if(0 <= x && x <= 3)
		return(-pow(x-3,3)/9.);
	return(0);
}

long double Basis2(long double x)
{
	if(0 <= x && x < 2)
		return(-pow(x,2)*(13.*x-36.)/72.);
	else if(0 <= x && x < 3)
		return(23.*pow(x,3)/72.-2.5*pow(x,2)+6.*x-4.);
	else if(0 <= x && x <= 4)
		return(-pow((x-4.)/2.,3));
	return(0);
}

long double Basis3(long double x)
{
	if(0 <= x && x < 2)
		return(pow(x,3)/24.);
	else if(0 <= x && x < 3)
		return(-3.*pow(x/2.,3)+2.5*pow(x,2)-5.*x+10./3.);
	else if(0 <= x && x < 4)
		return(11.*pow(x,3)/24.-5*pow(x,2)+17.5*x-115./6.);
	else if(0 <= x && x <= 5)
		return(-pow((x-5.),3)/6.);
	return(0);
}

long double Basisn(long double x)
{
	if(2 <= x && x < 3)
		return(pow(x-2.,3)/6.);
	else if(2 <= x && x < 4)
		return(-pow(x,3)/2.+5*pow(x,2)-16.*x+50./3.);
	else if(2 <= x && x < 5)
		return(pow(x,3)/2.-7*pow(x,2)+32.*x-142./3.);
	else if(2 <= x && x <= 6)
		return(-pow((x-6.),3)/6.);
	return(0);
}

void Import(char* File)
{
	ifstream Input(File);

	Input >> Coupling_Fraction >> Vacuum_Coupling;

	for(int i = 0; i < 259; i++)
		for(int j = 0; j < 151; j++)
			Input >> ImGf[j][i] >> ImGvCf[j][i] >> ImGvLf[j][i] >> ImGvQf[j][i] >> ImGVf[j][i] >> ReGvCf[j][i] >> ReGvLf[j][i] >> ReGvQf[j][i] >> ReGVf[j][i];

	return;
}
