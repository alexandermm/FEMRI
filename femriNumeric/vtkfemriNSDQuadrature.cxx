/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriNSDQuadrature.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/03 17:00:30 $
  Version:   $Revision: 1.3 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "vtkfemriNSDQuadrature.h"
#include "vtkObjectFactory.h"
#include "vtkCellType.h"
//alexmbcm: change: added the following line
#include "vtkMath.h"

vtkStandardNewMacro(vtkfemriNSDQuadrature);
vtkCxxRevisionMacro(vtkfemriNSDQuadrature, "$Revision: 1.3 $");

vtkfemriNSDQuadrature::vtkfemriNSDQuadrature()
{
  this->QuadraturePoints = NULL;
  this->QuadratureWeights = NULL;

  this->Order = 1;
  this->PreviousOrder = 0;
}

vtkfemriNSDQuadrature::~vtkfemriNSDQuadrature()
{
  if (this->QuadraturePoints)
  {
    this->QuadraturePoints->Delete();
    this->QuadraturePoints = NULL;
  }

  if (this->QuadratureWeights)
  {
    this->QuadratureWeights->Delete();
    this->QuadratureWeights = NULL;
  }
}

void vtkfemriNSDQuadrature::Initialize(int quadRuleName1, int quadRuleName2)
{
  if (this->QuadraturePoints)
  {
    this->QuadraturePoints->Delete();
    this->QuadraturePoints = NULL;
  }
  this->QuadraturePoints = vtkDoubleArray::New();
 
  if (this->QuadratureWeights)
  {
    this->QuadratureWeights->Delete();
    this->QuadratureWeights = NULL;
  }
  this->QuadratureWeights = vtkDoubleArray::New();
 
  
	vtkfemriNSDQuadrature* quadRule1 = vtkfemriNSDQuadrature::New();
	quadRule1->SetOrder(this->Order);
	vtkfemriNSDQuadrature* quadRule2 = vtkfemriNSDQuadrature::New();
	quadRule2->SetOrder(this->Order);
	  
	switch (quadRuleName1)
	{
		case GAUSS_LEGENDRE:
			quadRule1->Initialize1DLegendre();
			break;
        case GAUSS_LAGUERRE:
			quadRule1->Initialize1DLaguerre();
			break;
		case GAUSS_HERMITE:
			quadRule1->Initialize1DHermite();
			break;
        case GAUSS_HALF_HERMITE:
			quadRule1->Initialize1DHalfHermite();
			break;
	} 
	switch (quadRuleName2)
	{
		case GAUSS_LEGENDRE:
			quadRule2->Initialize1DLegendre();
			break;
        case GAUSS_LAGUERRE:
			quadRule2->Initialize1DLaguerre();
			break;
		case GAUSS_HERMITE:
			quadRule2->Initialize1DHermite();
			break;
        case GAUSS_HALF_HERMITE:
		quadRule2->Initialize1DHalfHermite();
			break;
	} 

	this->TensorProduct(quadRule1, quadRule2);
	  
	quadRule1->Delete();
	quadRule2->Delete();
}

void vtkfemriNSDQuadrature::ScaleTo01()
{
  if (this->QuadraturePoints->GetNumberOfComponents() != 1)
  {
    vtkErrorMacro("Error: scaling assumes Dimensionality == 1.");
    return;
  }
 
  double point[1];
  double weight;
  int numberOfQuadraturePoints = this->GetNumberOfQuadraturePoints();
  int i;
  for (i=0; i<numberOfQuadraturePoints; i++)
  {
    this->QuadraturePoints->GetTuple(i,point);
    point[0] = 0.5 * (point[0] + 1.0);
    this->QuadraturePoints->SetTuple(i,point);

    weight = this->QuadratureWeights->GetValue(i);
    weight *= 0.5;
    this->QuadratureWeights->SetValue(i,weight);
  }
}

//alexmbcm: changes: added the following inline functions used in the file 
//vtkThinPlateSplineTransform.cxx:

//------------------------------------------------------------------------
// some dull matrix things

inline double** vtkNewMatrix(int rows, int cols) 
{
  double *matrix = new double[rows*cols];
  double **m = new double *[rows];
  for(int i = 0; i < rows; i++) 
    {
    m[i] = &matrix[i*cols];
    }
  return m;
}

//------------------------------------------------------------------------
inline void vtkDeleteMatrix(double **m) 
{
  delete [] *m;
  delete [] m;
}

//------------------------------------------------------------------------
inline void vtkZeroMatrix(double **m, int rows, int cols) 
{
  for(int i = 0; i < rows; i++) 
    {
    for(int j = 0; j < cols; j++) 
      {
      m[i][j] = 0.0;
      }
    }
}

//------------------------------------------------------------------------

//This function finds the quadrature points and weights for 
//Gauss-Legendre quadrature
//Interval [0,1]
//Weight function: 1
void vtkfemriNSDQuadrature::Initialize1DLegendre()
{
  if (this->QuadraturePoints)
  {
    this->QuadraturePoints->Delete();
    this->QuadraturePoints = NULL;
  }
  this->QuadraturePoints = vtkDoubleArray::New();
 
  if (this->QuadratureWeights)
  {
    this->QuadratureWeights->Delete();
    this->QuadratureWeights = NULL;
  }
  this->QuadratureWeights = vtkDoubleArray::New();

	//alexmbcm: added Golub's method to get an arbitary number of gauss points
	if (this->Order == 0 || this->Order ==1)
	{
        this->QuadraturePoints ->SetNumberOfTuples(1);
		this->QuadratureWeights->SetNumberOfTuples(1);
		double* points  = this->QuadraturePoints->GetPointer(0);
		double* weights = this->QuadratureWeights->GetPointer(0);
		points[0]  = 0.5;
		weights[0] = 1.0;
		return;
	}
	
	int numPs;  //Number of points
	int i;      //Variable used for for loops
	//Zeroth moment of weight function 
	//(the integral of the weight function on quad interval) 
	double m_0 = 2.0; 
	double lengthSquared; //Variable used to calculate weight
  
    //Calculate max number of points needed based on max quad order 
	//needed
	if(this->Order % 2 == 0)
		numPs = (this->Order + 2)/2;
	else
		numPs = (this->Order + 1)/2;
	
	this->QuadraturePoints ->SetNumberOfTuples(numPs);
	this->QuadratureWeights->SetNumberOfTuples(numPs);
	
	//Fill points and weights using the Golub–Welsch algorithm
	double* points = this->QuadraturePoints->GetPointer(0);
	double* weights = this->QuadratureWeights->GetPointer(0);
    
	//Make matrix to apply Golub-Welsch algorithm
	//using the JacobiN function in vtkMath.h
	//Using inline functions above current function declaration
	double **A = vtkNewMatrix(numPs,numPs);
	vtkZeroMatrix(A,numPs,numPs);
		
	//Include in A the square root of 2 term recurrence coefficients
	//In the lower and upper diagonals
	for(i = 0; i < (numPs - 1); i++) 
		A[i+1][i] = A[i][i+1] = (i+1.0)/sqrt(4.0*(i+1.0)*(i+1.0) - 1.0);
		
	//Find points using JacobiN function
	double *qPoints       = new double[numPs];
	double **eigenVectors = vtkNewMatrix(numPs,numPs);
	vtkMath::JacobiN(A, numPs, qPoints, eigenVectors);
		
	//Put found points in point array taking into account 
	//they are ordered from positive to negative 
	for(i = 0; i < numPs; i++)
		points[i]  = qPoints[numPs-1 - i];
		
	//Calculate weights (first component of normalized eVect ^2 times fisrt moment of W)
	for(i = 0; i < numPs; i++)
		weights[i] = m_0*eigenVectors[0][numPs-1 - i]*eigenVectors[0][numPs-1 - i];
	
	//Delete matrices when done using them 
	vtkDeleteMatrix(A);
	delete [] qPoints; 
	vtkDeleteMatrix(eigenVectors);

  //Scale points and weights from [-1,1] to [0,1]
  this->ScaleTo01();
}

//This function finds the quadrature points and weights for 
//Gauss-Laguerre quadrature
//Interval [0,inf]
//Weight function: e^{-x}
void vtkfemriNSDQuadrature::Initialize1DLaguerre()
{
  if (this->QuadraturePoints)
  {
    this->QuadraturePoints->Delete();
    this->QuadraturePoints = NULL;
  }
  this->QuadraturePoints = vtkDoubleArray::New();
 
  if (this->QuadratureWeights)
  {
    this->QuadratureWeights->Delete();
    this->QuadratureWeights = NULL;
  }
  this->QuadratureWeights = vtkDoubleArray::New();
 
	//alexmbcm: added Golub's method to get an arbitary number of gauss points
	if (this->Order == 0 || this->Order ==1)
	{
        this->QuadraturePoints ->SetNumberOfTuples(1);
		this->QuadratureWeights->SetNumberOfTuples(1);
		double* points  = this->QuadraturePoints->GetPointer(0);
		double* weights = this->QuadratureWeights->GetPointer(0);
		points[0]  = 1.0;
		weights[0] = 1.0;
		return;
	}
	
	int numPs;  //Number of points
	int i;      //Variable used for for loops
	//Zeroth moment of weight function 
	//(just the integral of the weight function on quad interval) 
	double m_0 = 1.0; 
	double lengthSquared; //Variable used to calculate weight
  
    //Calculate max number of points needed based on max quad order 
	//needed
	if(this->Order % 2 == 0)
		numPs = (this->Order + 2)/2;
	else
		numPs = (this->Order + 1)/2;
	
	this->QuadraturePoints ->SetNumberOfTuples(numPs);
	this->QuadratureWeights->SetNumberOfTuples(numPs);
	
	//Fill points and weights using the Golub–Welsch algorithm
	double* points = this->QuadraturePoints->GetPointer(0);
	double* weights = this->QuadratureWeights->GetPointer(0);
    
	//Make 0 matrix to apply Golub-Welsch algorithm
	//using the JacobiN function in vtkMath.h
	//Using inline functions above current function declaration
	double **A = vtkNewMatrix(numPs,numPs);
	vtkZeroMatrix(A,numPs,numPs);
		
	//Include in A the square root of 2 term recurrence coefficients
	//In the lower and upper diagonals
	for(i = 0; i < numPs; i++) 
	{
		A[i][i] = 2.0*i + 1.0;
			
		if(i < numPs-1)
			A[i+1][i] = A[i][i+1] = i + 1.0;
	}
		
	//Find points using JacobiN function
	double *qPoints       = new double[numPs];
	double **eigenVectors = vtkNewMatrix(numPs,numPs);
	vtkMath::JacobiN(A, numPs, qPoints, eigenVectors);
		
	//Put found points (from 0 to 1) in point array taking into account 
	//they are ordered from positive to negative 
	for(i = 0; i < numPs; i++)
		points[i]  = qPoints[numPs-1 - i];
		
	//Calculate weights (first component of normalized eVect ^2 times first moment of W)
	for(i = 0; i < numPs; i++) 
		weights[i] = m_0*eigenVectors[0][numPs-1 - i]*eigenVectors[0][numPs-1 - i];
		
	//Delete matrices when done using them 
	vtkDeleteMatrix(A);
	delete [] qPoints; 
	vtkDeleteMatrix(eigenVectors);
}

//This function finds the quadrature points and weights for 
//Gauss-Hermite quadrature
//Interval [-inf,inf]
//Weight function: e^{-x^2}
void vtkfemriNSDQuadrature::Initialize1DHermite()
{
  if (this->QuadraturePoints)
  {
    this->QuadraturePoints->Delete();
    this->QuadraturePoints = NULL;
  }
  this->QuadraturePoints = vtkDoubleArray::New();
 
  if (this->QuadratureWeights)
  {
    this->QuadratureWeights->Delete();
    this->QuadratureWeights = NULL;
  }
  this->QuadratureWeights = vtkDoubleArray::New();
 
	//alexmbcm: added Golub's method to get an arbitary number of gauss points
	if (this->Order == 0 || this->Order ==1)
	{
        this->QuadraturePoints ->SetNumberOfTuples(1);
		this->QuadratureWeights->SetNumberOfTuples(1);
		double* points  = this->QuadraturePoints->GetPointer(0);
		double* weights = this->QuadratureWeights->GetPointer(0);
		points[0]  = 0.0;
		weights[0] = 1.772453850905516027298167483341;
		return;
	}
	
	int numPs;  //Number of points
	int i;      //Variable used for for loops
	//Zeroth moment of weight function 
	//(just the integral of the weight function on quad interval) 
	double m_0 = 1.772453850905516027298167483341; 
	double lengthSquared; //Variable used to calculate weight
  
    //Calculate max number of points needed based on max quad order 
	//needed
	if(this->Order % 2 == 0)
		numPs = (this->Order + 2)/2;
	else
		numPs = (this->Order + 1)/2;
	
	this->QuadraturePoints ->SetNumberOfTuples(numPs);
	this->QuadratureWeights->SetNumberOfTuples(numPs);
	
	//Fill points and weights using the Golub–Welsch algorithm
	double* points = this->QuadraturePoints->GetPointer(0);
	double* weights = this->QuadratureWeights->GetPointer(0);
    
	//Make 0 matrix to apply Golub-Welsch algorithm
	//using the JacobiN function in vtkMath.h
	//Using inline functions above current function declaration
	double **A = vtkNewMatrix(numPs,numPs);
	vtkZeroMatrix(A,numPs,numPs);
		
	//Include in A the square root of 2 term recurrence coefficients
	//In the lower and upper diagonals
	for(i = 0; i < (numPs - 1); i++) 
		A[i+1][i] = A[i][i+1] = sqrt((i+1.0) / 2.0);
		
	//Find points using JacobiN function
	double *qPoints       = new double[numPs];
	double **eigenVectors = vtkNewMatrix(numPs,numPs);
	vtkMath::JacobiN(A, numPs, qPoints, eigenVectors);
		
	//Put found points (from 0 to 1) in point array taking into account 
	//they are ordered from positive to negative 
	for(i = 0; i < numPs; i++)
		points[i]  = qPoints[numPs-1 - i];
		
	//Calculate weights (first component of normalized eVect ^2 times first moment of W)
	for(i = 0; i < numPs; i++) 
		weights[i] = m_0*eigenVectors[0][numPs-1 - i]*eigenVectors[0][numPs-1 - i];
		
	//Delete matrices when done using them 
	vtkDeleteMatrix(A);
	delete [] qPoints; 
	vtkDeleteMatrix(eigenVectors);
}

//This function finds the quadrature points and weights for 
//Gauss- half range Hermite quadrature
//Interval [0,inf]
//Weight function: e^{-x^2}
void vtkfemriNSDQuadrature::Initialize1DHalfHermite()
{
  if (this->QuadraturePoints)
  {
    this->QuadraturePoints->Delete();
    this->QuadraturePoints = NULL;
  }
  this->QuadraturePoints = vtkDoubleArray::New();
 
  if (this->QuadratureWeights)
  {
    this->QuadratureWeights->Delete();
    this->QuadratureWeights = NULL;
  }
  this->QuadratureWeights = vtkDoubleArray::New();
 
	//alexmbcm: added Golub's method to get an arbitary number of gauss points
	if (this->Order == 0 || this->Order ==1)
	{
        this->QuadraturePoints ->SetNumberOfTuples(1);
		this->QuadratureWeights->SetNumberOfTuples(1);
		double* points  = this->QuadraturePoints->GetPointer(0);
		double* weights = this->QuadratureWeights->GetPointer(0);
		points[0]  = 0.564189583547756286948079451561; //first alpha rec. coeff.
		weights[0] = 0.886226925452758013649083741671; //= m_0, see info below
		return;
	}

	int numPs;  //Number of points
	int i;      //Variable used for for loops
	//Zeroth moment of weight function 
	//(just the integral of the weight function on quad interval)
	//which is sqrt(pi)/2 
	double m_0 = 0.886226925452758013649083741671; 
	double lengthSquared; //Variable used to calculate weight
  
    //Calculate max number of points needed based on max quad order 
	//needed
	if(this->Order % 2 == 0)
		numPs = (this->Order + 2)/2;
	else
		numPs = (this->Order + 1)/2;
	
	//Check this since alpha and beta arrays are only 100 entries long
	if (numPs > 100)
	{
		vtkErrorMacro("femri Error: cannot generate half-range Hermite quadrature with >100 points");
		return;
	}

	
	this->QuadraturePoints ->SetNumberOfTuples(numPs);
	this->QuadratureWeights->SetNumberOfTuples(numPs);
	
	//Fill points and weights using the Golub–Welsch algorithm
	double* points = this->QuadraturePoints->GetPointer(0);
	double* weights = this->QuadratureWeights->GetPointer(0);
    
	//Make 0 matrix to apply Golub-Welsch algorithm
	//using the JacobiN function in vtkMath.h
	//Using inline functions above current function declaration
	double **A = vtkNewMatrix(numPs,numPs);
	vtkZeroMatrix(A,numPs,numPs);
		
	//Include in A the square root of 2 term recurrence coefficients
	//In the lower and upper diagonals
	double alpha[] =   {0.56418958354775628694807945156077,
						0.98842539284680028548706335878879,
						1.2859676193639399602827887260072,
						1.5247208440801153035130022763796,
						1.7301922743094392567715613980002,
						1.9134998431431025707186744531146,
						2.0806203364008332248176222241316,
						2.2352283805046391496583172950811,
						2.3797824435046374209405350458580,
						2.5160256434438664097634179698490,
						2.6452479250569531803261722004195,
						2.7684359535042559069132735915048,
						2.8863645940326945692708674993961,
						2.9996556533536035386906268743581,
						3.1088171759249201516926963722256,
						3.2142706360711282274489138373496,
						3.3163702970830873659208156975134,
						3.4154173324133389445368691110637,
						3.5116703446156295154051050586113,
						3.6053533459055664302954823302812,
						3.6966619115045907999685939739427,
						3.7857679927002249484946054222107,
						3.8728237301852214895392359107382,
						3.9579645104229992917022032498795,
						4.0413114410343916384468670119974,
						4.1229733747796281836030272921763,
						4.2030485788720019526602771862272,
						4.2816261227682035852100860580406,
						4.3587870403898888525030682138346,
						4.4346053100412971727359268038480,
						4.5091486858077933332742999188367,
						4.5824794070596267694605327867341,
						4.6546548072209951034250254583318,
						4.7257278387550172210545015603835,
						4.7957475280435453272958093455628,
						4.8647593712768612584458606016480,
						4.9328056804434948602898732480076,
						4.9999258868996881314362460126612,
						5.0661568087079578163326913880430,
						5.1315328868942965193196920906615,
						5.1960863949301972507737166643170,
						5.2598476250578105118341491173401,
						5.3228450545124346836883132956958,
						5.3851054942315484729257975795154,
						5.4466542222544300256576569466563,
						5.5075151036958845509104552954029,
						5.5677106989096931369032620838451,
						5.6272623612325201939687681911119,
						5.6861903255095047233332435598303,
						5.7445137884424375199090118491353,
						5.8022509816652991735273271228833,
						5.8594192383359449960986511943025,
						5.9160350539335612124739376867714,
						5.9721141418664674794919375700287,
						6.0276714844216699372382702394761,
						6.0827213795244305387450313134529,
						6.1372774837214812499446571090712,
						6.1913528517541013387434479172454,
						6.2449599730460277331169308290854,
						6.2981108053951899999360283032939,
						6.3508168061268025652438616206245,
						6.4030889609377744893289825207293,
						6.4549378106381760035349698739598,
						6.5063734759741768166952328786692,
						6.5574056806980573681919404101594,
						6.6080437730342609602385258146399,
						6.6582967456757199482476377461748,
						6.7081732544316109877349360217933,
						6.7576816356360647583633961233475,
						6.8068299224169953976698473146944,
						6.8556258599149692915076111589335,
						6.9040769195337678959671121791044,
						6.9521903122968986693306825350591,
						6.9999730013776709235245168072222,
						7.0474317138644914379712849613166,
						7.0945729518176711838144453628643,
						7.1414030026692022933628766113396,
						7.1879279490126046097749318465201,
						7.2341536778260020937955690657063,
						7.2800858891680256128175843237907,
						7.3257301043829101728031143857129,
						7.3710916738482261664035633942971,
						7.4161757842960244969119013675129,
						7.4609874657357568774760884437432,
						7.5055315980051307127752929190600,
						7.5498129169730510065279124767263,
						7.5938360204169703727139422368295,
						7.6376053735952952489577729747199,
						7.6811253145339664690624599020084,
						7.7244000590449317437148826745531,
						7.7674337054929440705959934148864,
						7.8102302393259426774548602123223,
						7.8527935373831919603671519420526,
						7.8951273719943601882384030760346,
						7.9372354148818055721842934739708,
						7.9791212408774955088030244620954,
						8.0207883314652089718801001570274,
						8.0622400781579563514236291414490,
						8.1034797857198902893237246604890,
						8.1445106752413705081427395333729 };

	double beta[] =    {0.88622692545275801364908374167057,
						0.18169011381620932846223247325497,
						0.34132512895943919856417178056476,
						0.50496215298800163193575115541755,
						0.67026419463961908567850839109470,
						0.83617049928031101554882352780077,
						1.0023478510110108422245382000471,
						1.1686711647442727438147851444475,
						1.3350829222423353579798779421644,
						1.5015525993447618438952914329038,
						1.6680623621881161688454808264457,
						1.8346010527937676419913598361001,
						2.0011613185512137843317104169486,
						2.1677381117632644853481959563471,
						2.3343278495405013980184823375142,
						2.5009279171337026699543207674084,
						2.6675363609572020882810741328780,
						2.8341516916678327579204958902571,
						3.0007727537827190275850360767873,
						3.1673986369644268117575670619156,
						3.3340286142031102452513292626711,
						3.5006620978281146517242161825195,
						3.6672986076183948893640172548364,
						3.8339377472958318506146765426260,
						4.0005791869361956805814791529341,
						4.1672226496283331945572198098899,
						4.3338679012299504436044302625654,
						4.5005147424120943373511796426753,
						4.6671630024168257031595689610456,
						4.8338125341123277420267271402514,
						5.0004632100412028349661109390409,
						5.1671149192366474464102065450405,
						5.3337675646378040182685821144089,
						5.5004210609766768917210652224104,
						5.6670753330391570642955682952632,
						5.8337303142250673602328783217237,
						6.0003859453488901957077630057946,
						6.1670421736354994723181628945983,
						6.3336989518748675754822847768167,
						6.5003562377071329380351548164824,
						6.6670139930151540677164102978028,
						6.8336721834061521384552925123128,
						7.0003307777675582697995981749345,
						7.1669897478849579708010267128251,
						7.3336490681122320857853129922995,
						7.5003087150857578832462612994777,
						7.6669686674759521591020735122514,
						7.8336289057705842156034430937754,
						8.0002894120852171919567543628648,
						8.1669501699968955260650680124847,
						8.3336111643978186628474948828657,
						8.5002723813662534111950360687225,
						8.6669338080523607985915613833777,
						8.8335954325769646773029096714318,
						9.0002572439415820669654687206362,
						9.1669192319482799641543432714415,
						9.3335813871281286949208497378615,
						9.5002437006771947561745369247992,
						9.6669061643991620951124682060914,
						9.8335687706537944739590941291968,
						10.000231512310556683634632403232,
						10.166894382706801955560032969015,
						10.333557375610009484376454465120,
						10.500220485183621585911135530124,
						10.666883705956086387928243278247,
						10.833547032792760509289109880864,
						11.000210460870368111787334097479,
						11.166873985653748998858438735472,
						11.333537602874659916843233494204,
						11.500201308512420585226118851013,
						11.666865098776219830209184429643,
						11.833528970088918017852670452953,
						12.000192919072200199562956622605,
						12.166856942532950351215969268374,
						12.333521037450731111772688580898,
						12.500185200966265767328353933291,
						12.666849430370830103998833752980,
						12.833513723096471358019471122426,
						13.000178076706979987205427563041,
						13.166842488889547515030463292031,
						13.333506957447050378213557768014,
						13.500171480290905645576587790805,
						13.666836055434449760707887593362,
						13.833500680986796172266830497575,
						14.000165355147131921918557887673,
						14.166830076199417020370969232692,
						14.333494842507453808659734643391,
						14.500159652510296519990481636418,
						14.666824504717973966732622886481,
						14.833489397707500712296539722582,
						15.000154330119154279097016194833,
						15.166819300652997918405552570938,
						15.333484308065630249211063224480,
						15.500149351167144682039038861734,
						15.666814428818282998384752156449,
						15.833479539927768773230251784075,
						16.000144683449807521410110924465,
						16.166809858381741531108640746823,
						16.333475063761848330821938603423,
						16.500140298667272629754611153557 };
  
	
	for(i = 0; i < numPs; i++) 
	{
		A[i][i] = alpha[i];
			
		if(i < numPs-1)
			A[i+1][i] = A[i][i+1] = sqrt(beta[i+1]);
	}
		
	//Find points using JacobiN function
	double *qPoints       = new double[numPs];
	double **eigenVectors = vtkNewMatrix(numPs,numPs);
	vtkMath::JacobiN(A, numPs, qPoints, eigenVectors);
		
	//Put found points (from 0 to 1) in point array taking into account 
	//they are ordered from positive to negative 
	for(i = 0; i < numPs; i++)
		points[i]  = qPoints[numPs-1 - i];
		
	//Calculate weights (first component of normalized eVect ^2 times first moment of W)
	for(i = 0; i < numPs; i++) 
		weights[i] = m_0*eigenVectors[0][numPs-1 - i]*eigenVectors[0][numPs-1 - i];
		
	//Delete matrices when done using them 
	vtkDeleteMatrix(A);
	delete [] qPoints; 
	vtkDeleteMatrix(eigenVectors);
}


//Make 2D product rule based on gauss quadratures used in NSD
void vtkfemriNSDQuadrature::TensorProduct(vtkfemriNSDQuadrature* quadRule1, vtkfemriNSDQuadrature* quadRule2)
{
  if (quadRule1->GetNumberOfQuadraturePoints() != quadRule2->GetNumberOfQuadraturePoints())
  {
    vtkErrorMacro("Error: cannot build tensor product rule if rules have different order.");
  }

  int numberOf1DQuadraturePoints = quadRule1->GetNumberOfQuadraturePoints();

  this->QuadraturePoints->SetNumberOfComponents(2);
  this->QuadraturePoints->SetNumberOfTuples(numberOf1DQuadraturePoints * numberOf1DQuadraturePoints);
  this->QuadratureWeights->SetNumberOfTuples(numberOf1DQuadraturePoints * numberOf1DQuadraturePoints);

  double point[2];
  double weight;
  int id;
  int i, j;
  
  // This loop computes the so called conical product rule for the standard triangle [0,0 1,0 0,1]
  // See also the code in LibMesh's src/quadrature/quadrature_conical.C
  for (j=0; j<numberOf1DQuadraturePoints; j++)
  {
    for (i=0; i<numberOf1DQuadraturePoints; i++)
    {
      id = j * numberOf1DQuadraturePoints + i;
        
      point[0] = quadRule1->GetQuadraturePoint(j)[0];
      point[1] = quadRule2->GetQuadraturePoint(i)[0];

      weight = quadRule1->GetQuadratureWeight(j) * quadRule2->GetQuadratureWeight(i);
      
      this->QuadraturePoints->SetTuple(id,point);
      this->QuadratureWeights->SetValue(id,weight);
    }
  }
}