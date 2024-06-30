#include "Expr.h"

std::vector<Variable> VariableList;

int set_Variable(std::string name, double value)
{
	Variable v;

	v.name = name;
	v.value = value;
	VariableList.push_back(v);
	return 0;
}

double expr_eval(const char* str, std::vector<double>& value)
{
	double x = 0.0;

	Expr exper(str);
	
	do{
		x = exper.expr();
		if (*exper.cur() == '\n') exper.next();
		value.push_back(x);
		if (exper.error() != EXPR_GOOD) break;
	} while (*exper.cur());

	if (exper.error())
	{
		printf("error [%s]\n", str);
		throw "error expression.";
	}
	return x;
}
double expr_eval(const char* str, std::vector<double*>& value)
{
	double x = 0.0;

	Expr exper(str);

	for (int i = 0; i < value.size(); i++)
	{
		x = exper.expr();
		double* p = value[i];
		*p = x;
		if (exper.error() != EXPR_GOOD) break;
	}
	if (*exper.cur() == '\n') exper.next();
	if (*exper.cur())
	{
		printf("error [%s] <-- garbage char[%c]\n", str, *exper.cur());
		exper.error(exper.cur(), 1, "garbage in expression");
		exper.error_chk();
	}
	if (exper.error())
	{
		printf("error [%s]\n", str);
		throw "error expression.";
	}

	return x;
}

double expr_eval(const char* str)
{


    double x = 0.0;

	Expr exper(str);

	if (*str == '%')
	{
		str++;
		char buf[80];
		char*p = buf;
		while (*str)
		{
			*p = *str;
			if (*str == '%') break;
			p++;
			str++;
		}
		if (*p == '%')
		{
			*p = '\0';
		}
		else
		{
			printf("error [%s]\n", str);
			throw "error expression( env value).";
		}

		double v;
		if (getenv(buf))
		{
			v = atof(getenv(buf));
			fprintf(stderr, "env =[%s]=%.16f\n", buf, v);
			printf("env =[%s]=%.16f\n", buf, v);
		}
		else
		{
			v = 0.0;
			fprintf(stderr, "endefined env =[%s]=NULL -> default value -> %.16f\n",  buf, v);
			printf( "endefined env =[%s]=NULL -> default value -> %.16f\n", buf, v);
		}
		return v;
	}
	x = exper.expr();
	if (*exper.cur() == '\n') exper.next();
	if (*exper.cur())
    {
		printf("error [%s] <-- garbage char[%c]\n", str, *exper.cur());
		exper.error(exper.cur(), 1, "garbage in expression");
		exper.error_chk();
	}
	if (exper.error())
	{
		printf("error [%s]\n", str);
		throw "error expression.";
	}
	return x;
}

