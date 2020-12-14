#include <cmath>

template <typename F>
struct Vec3{
	F x, y, z;

#if 1
	template <typename G>
	operator Vec3<G>() const {
		return {
			static_cast<G>(x),
			static_cast<G>(y),
			static_cast<G>(z)
		};
	}
#endif
#if 0
	template <typename G>
	Vec3 &operator=(const Vec3<G> &src) & {
		return *this = {
			static_cast<F>(x),
			static_cast<F>(y),
			static_cast<F>(z)
		};
	}
#endif
};


// a * b + c without intermediate rounding
template <typename F>
inline F tp_fma(F a, F b, F c);

template<>
inline float tp_fma(float a, float b, float c){ return ::fmaf(a, b, c); }

template<>
inline double tp_fma(double a, double b, double c){ return ::fma(a, b, c); }

template <typename F>
inline F prec_cross_prod_inner(F ux, F uy, F vx, F vy){
	F p = ux * vy;
	F r = tp_fma(ux, vy, -p);
	return r + tp_fma(-vx, uy, p);
}

template <typename F>
__attribute__((noinline))
Vec3<F> prec_cross_prod(Vec3<F> u, Vec3<F> v){
	return{
		prec_cross_prod_inner(u.y, u.z, v.y, v.z),
		prec_cross_prod_inner(u.z, u.x, v.z, v.x),
		prec_cross_prod_inner(u.x, u.y, v.x, v.y)
	};
}

template <typename F>
__attribute__((noinline))
Vec3<F> cross_prod(Vec3<F> u, Vec3<F> v){
	return{
		(u.y * v.z - u.z * v.y),
		(u.z * v.x - u.x * v.z),
		(u.x * v.y - u.y * v.x)
	};
}

#include <cstdio>
#include <cstdlib>

template <typename F>
int vec_print(
		const Vec3<F> &v,
		const char *form = "(%e, %e, %e)\n",
		FILE *fp = stdout)
{
	return fprintf(fp, form, v.x, v.y, v.z);
}

int main(int ac, char **av){
	srand48(ac>1 ? atof(av[1]) : 20201212);

	Vec3<double> u, v;

	auto rnd = []{
		return drand48() - 0.5;
	};

	u.x = rnd() + 0.5;
	u.y = rnd() + 0.5;
	u.z = rnd() + 0.5;

	const double eps = 1.0e-7;

	// 微小回転でほぼ平行なベクトルを作る
	v.x = u.x + eps * (rnd() * u.y + rnd() * u.z);
	v.y = u.y + eps * (rnd() * u.z + rnd() * u.x);
	v.z = u.z + eps * (rnd() * u.x + rnd() * u.y);

	Vec3<float>  uf = u,  vf = v;
	Vec3<double> ud = uf, vd = vf;

	auto wf0 = cross_prod      (uf, vf);
	auto wf  = prec_cross_prod (uf, vf);
	auto wd  = cross_prod      (ud, vd);
	auto wd1 = prec_cross_prod (ud, vd);
	auto wd2=  prec_cross_prod (u,  v ); // 入力を単精度に丸める前のもの

	const char *form1 =  "(%24.16e, %24.16e, %24.16e)\n";
	vec_print(wf0, form1);
	vec_print(wf,  form1);
	vec_print(wd,  form1);
	vec_print(wd1, form1);
	vec_print(wd2, form1);

	const char *form2 =  "(%A, %A, %A)\n";
	vec_print(wf0, form2);
	vec_print(wf,  form2);
	vec_print(wd,  form2);
	vec_print(wd1, form2);
	vec_print(wd2, form2);

	return 0;
}

