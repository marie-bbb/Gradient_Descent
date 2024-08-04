#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath> 
#include <sstream>
#include <chrono>
#include <vector>
#include <string>

using namespace std;

const int initialTime(0);
const int finalTime(1);
const double timeStep(0.05);
double a = 0.0;
const int maxGdIt = 3000;
const double delta = 0.00009;
const double alpha = 0.01;

int A = 0;
int B = 0;

/*signum function*/
double sign(double n)
{
	if (n > 0.0)
		return 1.0;
	else if (n < 0.0)
		return -1.0;
	else if (n == 0.0)
		return 0.0;
}

class GradientDescentMethod
{
public:
	GradientDescentMethod(int sizeX, int sizeY, double t, double dt);															
	int GetDegreesOfFreedom();
	void SetW(int& sizeX);																	
	void SetH(int& sizeY);																	
	void SetK(int k);																		
	void SetHX(double _hx);																	
	void SetHY(double _hy);
	void SetT(double t);
	void SetDT(double _dt);
	int GetW() const;																		
	int GetH() const;																		
	int GetK() const;																		
	double GetHX() const;
	double GetHY() const;
	double GetT() const;
	double GetDT() const;

	void PrintFunc1D(double* f, const char* filename, int gradDes);
	void PrintFunc2D(double* f1, double* f2, const char* filename, int gradDes);

	void SetInitGamma(double* gamma1, double* u, double* I, int gradDes);		
	void SetInitFunction(double* f, int gradDes, const char* function);						
	void SetInitData(double* f, double* data);												
	void SetInitA(double* a1, double* a2);										

	void PrimaryU(double* u, double* a1, double* a2, double* data, const int& t, const char* filename, int gradDes);	
	void AdjointGamma(double* u, double* gamma1, double* a1, double* a2, double* I, const int& t, int gradDes);								
	void VelocityUpdate(double* u, double* gamma1, double* a1, double* a2, const double& delta, const double& aplha, int gradDes);

	~GradientDescentMethod() {}

private:
	int W;			//image width
	int H;			//image height
	double hx;		//step on the x axis
	double hy;		//step on the y axis 

	double T, dt;	//real time, time step
	int K;			//Ñumber of images
};

GradientDescentMethod::GradientDescentMethod(int sizeX, int sizeY, double t, double dt)
{
	SetW(sizeX);
	SetH(sizeY);
	SetHX(1.0 / (double)W);
	SetHY(1.0 / (double)H);
	SetT(t);
	SetDT(2.0 * this->hx);
	SetK(ceil(T / dt));
}

void GradientDescentMethod::SetW(int& sizeX)
{
	this->W = sizeX;
}
void GradientDescentMethod::SetH(int& sizeY)
{
	this->H = sizeY;
}
void GradientDescentMethod::SetK(int k)
{
	this->K = k;
}
void GradientDescentMethod::SetHX(double _hx)
{
	this->hx = _hx;
}
void GradientDescentMethod::SetHY(double _hy)
{
	this->hy = _hy;
}
void GradientDescentMethod::SetT(double t)
{
	this->T = t;
}
void GradientDescentMethod::SetDT(double _dt)
{
	this->dt = _dt;
}
int GradientDescentMethod::GetW() const
{
	return this->W;
}
int GradientDescentMethod::GetH() const
{
	return this->H;
}
int GradientDescentMethod::GetK() const
{
	return this->K;
}
double GradientDescentMethod::GetHX() const
{
	return this->hx;
}
double GradientDescentMethod::GetHY() const
{
	return this->hy;
}
double GradientDescentMethod::GetT() const
{
	return this->T;
}
double GradientDescentMethod::GetDT() const
{
	return this->dt;
}

int GradientDescentMethod::GetDegreesOfFreedom()
{
	return this->W + this->H * this->W + this->W * this->H * (this->K + 1);
}

/*Save the scalar data*/
void GradientDescentMethod::PrintFunc1D(double* f, const char* filename, int gradDes)
{
	for (int k = 1; k <= K; k++)
	{
		if ((k % 10 == 0 && (gradDes % 10 == 0 || (gradDes > 860 && gradDes < 870))) || filename == "kolecko")
		{
			stringstream s;
			s << filename << "-" << k << "_" << gradDes << ".txt";
			fstream file1;

			file1.open(s.str(), ios::out);
			for (int j = 1; j <= H; j++)
			{
				for (int i = 1; i <= W; i++)
				{
					int p = i + j * W + k * W * H;
					file1 << i * hx << "\t" << j * hy << "\t \t" << f[p] << endl;
				}
				file1 << endl;
			}
			file1.close();
		}
	}
}

/*Save the vector data*/
void GradientDescentMethod::PrintFunc2D(double* f1, double* f2, const char* filename, int gradDes)
{
	for (int k = 1; k <= K; k++)
	{
		if ((k % 10 == 0 && (gradDes % 10 == 0 || (gradDes > 860 && gradDes < 870))) || filename == "kolecko")
		{
			stringstream s;
			s << filename << "-" << k << "_" << gradDes << ".txt";
			fstream file1;

			file1.open(s.str(), ios::out);
			for (int j = 1; j <= H; j++)
			{
				for (int i = 1; i <= W; i++)
				{
					int p = i + j * W + k * W * H;
					file1 << i * hx << "\t" << j * hy << "\t" << f1[p] << "\t \t" << f2[p] << endl;
				}
				file1 << endl;
			}
			file1.close();
		}
	}
}

/*Get dimensions of the files*/
int ReadDimPGM(const char* fileName)
{
	std::fstream file;
	file.open(fileName, std::ios::in);
	if (!file)
	{
		std::cerr << "Unable to open the file " << fileName << std::endl;
		return false;
	}

	std::string magicNumber;
	file >> magicNumber;
	if (magicNumber != "P2" && magicNumber != "P5")
	{
		std::cerr << "Unsupported format in file " << fileName << ". Only PGM format is supported." << std::endl;
		return false;
	}
	bool binary = (magicNumber == "P5");

	char character;
	file.get(character);
	while (!file.eof() and (character == ' ' || character == '\t' || character == '\r' || character == '\n'))
	{
		file.get(character);
		if (character == '#')
			while (!file.eof() && (character != '\n'))
				file.get(character);
	}
	file.unget();

	int maxColors;
	file >> A >> B >> maxColors;

	return A * B * A * B;
}

/*Read the picture files*/
bool ReadPGM(GradientDescentMethod problem, const char* fileName, int width, double* data, int m)
{
	std::fstream file;
	file.open(fileName, std::ios::in);
	if (!file)
	{
		std::cerr << "Unable to open the file " << fileName << std::endl;
		return false;
	}

	std::string magicNumber;
	file >> magicNumber;
	if (magicNumber != "P2" && magicNumber != "P5")
	{
		std::cerr << "Unsupported format in file " << fileName << ". Only PGM format is supported." << std::endl;
		return false;
	}
	bool binary = (magicNumber == "P5");

	char character;
	file.get(character);
	while (!file.eof() and (character == ' ' || character == '\t' || character == '\r' || character == '\n'))
	{
		file.get(character);
		if (character == '#')
			while (!file.eof() && (character != '\n'))
				file.get(character);
	}
	file.unget();

	int maxColors, w, h;
	file >> w >> h >> maxColors;

	for (int j = 1; j <= h; j++)
	{
		for (int i = 1; i <= w; i++)
		{
			int col;
			unsigned char col_aux;
			if (binary)
			{
				file >> col_aux;
				col = (int)col_aux;
			}
			else file >> col;
			int s = i + j * problem.GetW() + m * problem.GetW() * problem.GetH();
			data[s] = (double)col / (double)maxColors;
		}
	}
	problem.PrintFunc1D(data, "kolecko", 0);
	return true;
}

/*Initial condition - Gauss or signum function*/
void GradientDescentMethod::SetInitFunction(double* f, int gradDes, const char* function)
{
	for (int k = 1; k <= K; k++)
	{
		for (int j = 1; j <= H; j++)
		{
			for (int i = 1; i <= W; i++)
			{
				double x = i * hx;
				double y = j * hy;
				if (function == "gauss")
				{
					f[i + j * W + k * W * H] = exp((-((x - 10.0) * (x - 10.0)) - ((y - 10.0) * (y - 10.0))) / 2.0);
				}
				else if (function == "signum")
				{
					double v = (x - 10.0) * (x - 10.0) + (y - 10.0) * (y - 10.0) - 10.0;
					f[i + j * W + k * W * H] = sign(v) / 2.0 - 0.5;
				}
				else
				{
					cout << "Function not available for initialization." << endl;
				}
			}
		}
	}
}

/*Initial condition*/
void GradientDescentMethod::SetInitData(double* f, double* data)
{
	for (int k = 1; k <= K; k++)
	{
		for (int j = 1; j <= H; j++)
		{
			for (int i = 1; i <= W; i++)
			{
				f[i + j * W + k * W * H] = data[i + j * W + k * W * H];
			}
		}
	}
}

/*Initial condition - adjoint problem*/
void GradientDescentMethod::SetInitGamma(double* gamma1, double* u, double* I, int gradDes)
{
	for (int k = 1; k <= K; k++)
	{
		for (int j = 1; j <= H; j++)
		{
			for (int i = 1; i <= W; i++)
			{
				int r = i + j * W + k * W * H;
				gamma1[r] = 0.0;
			}
		}
	}
	PrintFunc1D(gamma1, "outputg0", gradDes);
}

/*Initial estimate of optical flow*/
void GradientDescentMethod::SetInitA(double* a1, double* a2)
{
	for (int k = 1; k <= K; k++)
	{
		for (int j = 1; j <= H; j++)
		{
			for (int i = 1; i <= W; i++)
			{
				a1[i + j * W + k * W * H] = a;
				a2[i + j * W + k * W * H] = a;
			}
		}
	}
}


/*Primary problem*/
void GradientDescentMethod::PrimaryU(double* u, double* a1, double* a2, double* data, const int& t, const char* filename, int gradDes)
{
	SetInitData(u, data);

	for (int k = 1; k <= K; k++)
	{
		for (int j = 1; j <= H; j++)
		{
			for (int i = 1; i <= W; i++)
			{
				int s = i + j * W + k * W * H;

				if (i == 1 || j == 1 || i == (W) || j == (H) || k == 1)
				{
					u[s + W * H] = u[s];
				}
				else
				{
					if (a1[s] < 0 && a2[s] < 0)
					{
						u[s + W * H] = u[s] - ((a1[s] * dt) / hx) * (u[s + 1] - u[s]) - ((a2[s] * dt) / hy) * (u[s + W] - u[s]);
					}
					else if (a1[s] >= 0 && a2[s] < 0)
					{
						u[s + W * H] = u[s] - ((a1[s] * dt) / hx) * (u[s] - u[s - 1]) - ((a2[s] * dt) / hy) * (u[s + W] - u[s]);
					}
					else if (a1[s] < 0 && a2[s] >= 0)
					{
						u[s + W * H] = u[s] - ((a1[s] * dt) / hx) * (u[s + 1] - u[s]) - ((a2[s] * dt) / hy) * (u[s] - u[s - W]);
					}
					else
					{
						u[s + W * H] = u[s] - ((a1[s] * dt) / hx) * (u[s] - u[s - 1]) - ((a2[s] * dt) / hy) * (u[s] - u[s - W]);
					}
				}
			}
		}
	}
	PrintFunc1D(u, filename, gradDes);
}

/*Adjoint problem*/
void GradientDescentMethod::AdjointGamma(double* u, double* gamma1, double* a1, double* a2, double* I, const int& t, int gradDes)
{
	SetInitGamma(gamma1, u, I, gradDes);

	for (int k = 1; k <= K; k++)
	{
		for (int j = 1; j <= H; j++)
		{
			for (int i = 1; i <= W; i++)
			{
				int r = i + j * W + k * W * H;

				if (i == 1 || i == (W) || j == 1 || j == (H))
				{
					gamma1[r + W * H] = 0.0;
				}
				else
				{
					if (a1[r] < 0.0 && a2[r] < 0.0)
					{
						gamma1[r + W * H] = gamma1[r] - ((a1[r] * dt) / hx) * (gamma1[r + 1] - gamma1[r]) - ((a2[r] * dt) / hy) * (gamma1[r + W] - gamma1[r])
							+ dt * 2.0 * (I[r] - u[r]);
					}
					else if (a1[r] >= 0.0 && a2[r] < 0.0)
					{
						gamma1[r + W * H] = gamma1[r] - ((a1[r] * dt) / hx) * (gamma1[r] - gamma1[r - 1]) - ((a2[r] * dt) / hy) * (gamma1[r + W] - gamma1[r])
							+ dt * 2.0 * (I[r] - u[r]);
					}
					else if (a1[r] < 0.0 && a2[r] >= 0.0)
					{
						gamma1[r + W * H] = gamma1[r] - ((a1[r] * dt) / hx) * (gamma1[r + 1] - gamma1[r]) - ((a2[r] * dt) / hy) * (gamma1[r] - gamma1[r - W])
							+ dt * 2.0 * (I[r] - u[r]);
					}
					else
					{
						gamma1[r + W * H] = gamma1[r] - ((a1[r] * dt) / hx) * (gamma1[r] - gamma1[r - 1]) - ((a2[r] * dt) / hy) * (gamma1[r] - gamma1[r - W])
							+ dt * 2.0 * (I[r] - u[r]);
					}
				}
			}
		}
	}
	PrintFunc1D(gamma1, "outputg", gradDes);
}

/*Updating the optical flow - gradient descent*/
void GradientDescentMethod::VelocityUpdate(double* u, double* gamma1, double* a1, double* a2, const double& delta, const double& alph, int gradDes)
{
	for (int k = 1; k <= K; k++)
	{
		for (int j = 1; j <= H; j++)
		{
			for (int i = 1; i <= W; i++)
			{
				int p = i + j * W + k * W * H;
				a1[p] = a1[p] - delta * gamma1[p] * ((u[p + 1] - u[p]) / hx); //- alph * a1[p];
				a2[p] = a2[p] - delta * gamma1[p] * ((u[p + W] - u[p]) / hy); //- alph * a2[p];
			}
		}
	}
	PrintFunc2D(a1, a2, "outputa", gradDes);
}


int main()
{
	system("CHCP 1250 > NUL");

	vector<const char*> filenames =
	{
		"circle1.pgm", "circle2.pgm", "circle3.pgm", "circle4.pgm",	"circle5.pgm", 
		"circle6.pgm", "circle7.pgm", "circle8.pgm", "circle9.pgm", "circle10.pgm",
		"circle11.pgm", "circle12.pgm", "circle13.pgm", "circle14.pgm", "circle15.pgm",
		"circle16.pgm",	"circle17.pgm", "circle18.pgm", "circle19.pgm", "circle20.pgm"
	};

	vector<ifstream> files(filenames.size());

	for (size_t k = 0; k < filenames.size(); k++)
	{
		files[k].open(filenames[k]);
	}

	int C = ReadDimPGM(filenames[0]);
	int D = C / (A * B * B);
	int E = C / (B * A * A);
	cout << D << " " << E << endl;

	GradientDescentMethod problem(D, E, finalTime, timeStep);

	double* u = new double[problem.GetDegreesOfFreedom()];
	double* gamma1 = new double[problem.GetDegreesOfFreedom()];
	double* a1 = new double[problem.GetDegreesOfFreedom()];
	double* a2 = new double[problem.GetDegreesOfFreedom()];
	double* I = new double[problem.GetDegreesOfFreedom()];

	for (int k = 0; k < 20; k++)
	{
		ReadPGM(problem, filenames[k], C, I, k + 1);
	}

	/*Initialize the velocity vector*/
	problem.SetInitA(a1, a2);

	/*Gradient descent iteration*/
	for (int gdIt = 1; gdIt <= maxGdIt; gdIt++)
	{
		problem.PrimaryU(u, a1, a2, I, finalTime, "outputu", gdIt);
		problem.AdjointGamma(u, gamma1, a1, a2, I, finalTime, gdIt);
		problem.VelocityUpdate(u, gamma1, a1, a2, delta, alpha, gdIt);

		if (gdIt % 10 == 0)
		{
			cout << "Konec iterace èíslo " << gdIt << endl;
			cout << "*______________________*" << endl;
		}
	}
}