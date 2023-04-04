#ifndef PARTICLELIST_H
#define PARTICLELIST_H

#include <vector>
#include "particle.hpp"


// List of	confirmed resonances (from the PDG 2022 edition) made of up, down,
// and strange quarks.  Resonances with charm and heavier are excluded.



/*-----------------------------------------------------------------------------
 Global data used in this program.  Listed are confirmed resonances from
 2022 pdg.  (The mesons are listed in pdg meson summary starting pg 34.
 A concise list of mesons used are mesons on pg 77 in the table with a dot by
 them--indicating existence is confirmed.  The baryons used are those from the
 table on page 78 with 3 or 4 stars--indicating confirmation.  Baryons with 2
 or fewer stars were not used.  Baryon properties were taken from page 79 and
 on of pdg baryon summary.) This selection strategy was used in
 arXiv:0901.1430v1.

 When masses are listed as a range, like 1200-1500, the midpoint is used.
 When degeneracies are not listed, I make up conservative (ie, small) values
 so I underestimate instead of overestimate the resonance contribution.


 Entries below are of the form:
 (name, e mass, J_degen, Fs, statistics)
 J_degen = 2J+1,
 Fs is a vector of flavor states which are given a triplet of
 quantum numbers {Baryon Number, Charge, Strangeness (,Charm)}
*/

//Entries with //! are new since the 2012 listing

const std::vector<Particle> ParticleList = {

	//-----------------------------------------------------------------------------
	//Light Unflavored Mesons:

	Particle("pi+", 139.57, 1, {{0,1,0},{0,-1,0}}, BOSON),
	Particle("pi0", 134.98, 1, {{0,0,0}}, BOSON),
	Particle("eta", 547.86, 1, F_Iscl, BOSON),

	//Many sources omit this resonance due to its huge width
	//Particle("f0(500)", 475.0, 1, 1, 0, 0, BOSON),

	Particle("rho(770)", 775.3, 3, F_Ivect, BOSON),
	Particle("omega(782)", 782.7, 3, F_Iscl, BOSON),
	Particle("eta'(958)", 957.8, 1, F_Iscl, BOSON),
	Particle("f0(980)", 990.0, 1, F_Iscl, BOSON),
	Particle("a0(980)", 980.0, 1, F_Ivect, BOSON),
	Particle("phi(1020)", 1019.5, 3, F_Iscl, BOSON),

	Particle("h1(1170)", 1166.0, 3, F_Iscl, BOSON),
	Particle("b1(1235)", 1229.5, 3, F_Ivect, BOSON),
	Particle("a1(1260)", 1230.0, 3, F_Ivect, BOSON),
	Particle("f2(1270)", 1275.5, 5, F_Iscl, BOSON),
	Particle("f1(1285)", 1281.9, 3, F_Iscl, BOSON),
	Particle("eta(1295)", 1294.0, 1, F_Iscl, BOSON),
	Particle("pi(1300)", 1300.0, 1, F_Ivect, BOSON),
	Particle("a2(1310)", 1318.2, 5, F_Ivect, BOSON),
	Particle("f0(1370)", 1350.0, 1, F_Iscl, BOSON),
	Particle("pi1(1400)", 1354.0, 3, F_Ivect, BOSON),
	Particle("eta(1405)", 1408.8, 1, F_Iscl, BOSON),
	Particle("h1(1415)", 1416.0, 3, F_Iscl, BOSON), //!
	Particle("f1(1420)", 1426.3, 3, F_Iscl, BOSON),
	Particle("omega(1420)", 1410.0, 3, F_Iscl, BOSON),
	Particle("a0(1450)", 1474.0, 1, F_Ivect, BOSON),
	Particle("rho(1450)", 1465.0, 3, F_Ivect, BOSON),
	Particle("eta(1475)", 1475.0, 1, F_Iscl, BOSON),
	Particle("f0(1500)", 1506.0, 1, F_Iscl, BOSON),
	Particle("f2'(1525)", 1517.4, 5, F_Iscl, BOSON),
	Particle("pi1(1600)", 1661.0, 3, F_Ivect, BOSON),
	Particle("a1(1640)", 1655.0, 3, F_Ivect, BOSON), //!
	Particle("eta2(1645)", 1617.0, 5, F_Iscl, BOSON),
	Particle("omega(1650)", 1670.0, 3, F_Iscl, BOSON),
	Particle("omega3(1670)", 1667.0, 7, F_Iscl, BOSON),
	Particle("pi2(1670)", 1670.6, 5, F_Ivect, BOSON),
	Particle("phi(1680)", 1680, 3, F_Iscl, BOSON),
	Particle("rho3(1690)", 1688.8, 7, F_Ivect, BOSON),
	Particle("rho(1700)",  1720.0, 3, F_Ivect, BOSON),
	Particle("a2(1700)",  1698.0, 5, F_Ivect, BOSON), //!
	Particle("f0(1710)", 1704.0, 1, F_Iscl, BOSON),
	Particle("pi(1800)", 1810.0, 1, F_Ivect, BOSON),
	Particle("phi3(1850)", 1854.0, 7, F_Iscl, BOSON),
	Particle("eta2(1870)", 1842.0, 5, F_Iscl, BOSON), //!
	Particle("pi2(1880)", 1874.0, 5, F_Ivect, BOSON),
	Particle("f2(1950)" , 1936.0, 5, F_Iscl, BOSON),
	Particle("a4(1970)", 1967.0, 9, F_Ivect, BOSON),

	//Passing 2 GeV threshold--some papers end here
	Particle("f2(2010)", 2011.0, 5, F_Iscl, BOSON),
	Particle("f4(2050)", 2018.0, 9, F_Iscl, BOSON),
	Particle("phi(2170)", 2162.0, 3, F_Iscl, BOSON),
	Particle("f2(2300)", 2297.0, 5, F_Iscl, BOSON),
	Particle("f2(2340)", 2345.0, 5, F_Iscl, BOSON),


	//-----------------------------------------------------------------------------
	//Strange Mesons (pdg starting on pg 39):

	Particle("K+", 493.68, 1, {{0,1,1},{0,-1,-1}}, BOSON),
	Particle("K0", 497.61, 1, {{0,1,0},{0,-1,0}}, BOSON),

	//I do not include KL or KS since they are not independent of K0 and anti-K0

	Particle("K*(700)", 845.0, 1, F_K, BOSON), //!
	Particle("K*(892)+", 891.7, 3, {{0,1,1},{0,-1,-1}}, BOSON),
	Particle("K*(892)0", 895.6, 3, {{0,1,0},{0,-1,0}}, BOSON),
	Particle("K*2(1430)+", 1427.3, 5, {{0,1,1},{0,-1,-1}}, BOSON),
	Particle("K*2(1430)0", 1432.4, 5, {{0,1,0},{0,-1,0}}, BOSON),
	Particle("K1(1270)", 1253.0, 3, F_K, BOSON),
	Particle("K1(1400)", 1403.0, 3, F_K, BOSON),
	Particle("K*(1410)", 1414.0, 3, F_K, BOSON),
	Particle("K*0(1430)", 1425.0, 1, F_K, BOSON),
	Particle("K*(1680)", 1718.0, 3, F_K, BOSON),
	Particle("K2(1770)", 1773.0, 5, F_K, BOSON),
	Particle("K3*(1780)", 1779.0, 7, F_K, BOSON),
	Particle("K2(1820)", 1819.0, 5, F_K, BOSON),
	Particle("K2*(1980)", 1994.0, 5, F_K, BOSON), //!

	//Passing 2 GeV threshold--some papers end here
	Particle("K4*(2045)", 2048.0, 9, F_K, BOSON),


	/*  Exclude particles with charm:
	//-----------------------------------------------------------------------------
	//Charmed Mesons (pdg starting on pg 43):

	Particle("D+", 1869.7, 1, F_D1, BOSON),
	Particle("D0", 1864.8, 1, F_D0, BOSON),

	//CHECK THESE:
	Particle("D*(2007)0", 2006.9, 3, F_D0, BOSON),
	Particle("D*(2010)+", 2010.3, 3, F_D1, BOSON),
	Particle("D*(2300)0", 2343.0, 1, F_D0, BOSON),

	// D*(2300)+  not confirmed, not listed here //!

	Particle("D1(2420)0", 2420.1, 3, F_D0, BOSON),
	Particle("D1(2420)+", 2424.1, 3, F_D1, BOSON),
	Particle("D1(2430)0", 2412.0, 3, F_D0, BOSON), //!

	Particle("D2*(2460)0", 2459.9, 5, F_D0, BOSON),
	Particle("D2*(2460)+", 2462.3, 5, F_D1, BOSON),

	//Definitely exclude these
	Particle("D3*(2750)0", 2763.1, 7, F_D0, BOSON), //!
	Particle("D3*(2750)+", 2763.1, 7, F_D1, BOSON), //!

	//-----------------------------------------------------------------------------
	//Charmed, Strange Mesons (pdg starting on pg 48):

	//careful with the degeneracy: must double pdg value to include anti-particle

	Particle("Ds+", 1968.4, 1, F_Ds, BOSON),

	//This is not certain: J is not listed, so I use value J=1
	Particle("Ds*+", 2112.2, 3, F_Ds, BOSON),

	Particle("Ds0*(2317)+", 2317.8, 1, F_Ds, BOSON),
	Particle("Ds1*(2460)+", 2459.5, 3, F_Ds, BOSON),
	Particle("Ds1(2536)+", 2535.1, 3, F_Ds, BOSON),
	Particle("Ds2*(2573)+", 2569.1, 5, F_Ds, BOSON),

	Particle("Ds1*(2700)+", 2714.1, 3, F_Ds, BOSON), //!
	Particle("Ds3*(2860)+", 2860.0, 7, F_Ds, BOSON), //!

	//-----------------------------------------------------------------------------
	//c-cbar Mesons:

	Particle("etac(1S)", 2983.9, 1, F_Iscl, BOSON),

	*/
	//-----------------------------------------------------------------------------
	//N Baryons (pgd starting on page 79):

	Particle("p", 938.27, 2, {{1,1,0},{-1,-1,0}}, FERMION),
	Particle("n", 939.57, 2, {{1,0,0},{-1,0,0}}, FERMION),

	Particle("N(1440)", 1440.0, 2, F_N, FERMION),
	Particle("N(1520)", 1515.0, 4, F_N, FERMION),
	Particle("N(1535)", 1530.0, 2, F_N, FERMION),
	Particle("N(1650)", 1650.0, 2, F_N, FERMION),
	Particle("N(1675)", 1675.0, 6, F_N, FERMION),
	Particle("N(1680)", 1685.0, 6, F_N, FERMION),
	Particle("N(1700)", 1720.0, 4, F_N, FERMION),
	Particle("N(1710)", 1710.0, 2, F_N, FERMION),
	Particle("N(1720)", 1720.0, 4, F_N, FERMION),
	Particle("N(1875)", 1875.0, 4, F_N, FERMION),
	Particle("N(1880)", 1880.0, 2, F_N, FERMION), //!
	Particle("N(1895)", 1895.0, 2, F_N, FERMION), //!
	Particle("N(1900)", 1920.0, 4, F_N, FERMION),

	Particle("N(2060)", 2100.0, 6, F_N, FERMION), //!
	Particle("N(2100)", 2100.0, 2, F_N, FERMION), //!
	Particle("N(2120)", 2120.0, 4, F_N, FERMION), //!
	Particle("N(2190)", 2180.0, 8, F_N, FERMION),
	Particle("N(2220)", 2250.0, 10, F_N, FERMION),
	Particle("N(2250)", 2280.0, 10, F_N, FERMION),
	Particle("N(2600)", 2600.0, 12, F_N, FERMION),


	//-----------------------------------------------------------------------------
	//Delta Baryons (pgd starting on pg 82)

	Particle("Delta(1232)", 1232.0, 4, F_Del, FERMION),
	Particle("Delta(1600)", 1570.0, 4, F_Del, FERMION),
	Particle("Delta(1620)", 1610.0, 2, F_Del, FERMION),
	Particle("Delta(1700)", 1710.0, 4, F_Del, FERMION),
	Particle("Delta(1900)", 1860.0, 2, F_Del, FERMION), //!
	Particle("Delta(1905)", 1880.0, 6, F_Del, FERMION),
	Particle("Delta(1910)", 1900.0, 2, F_Del, FERMION),
	Particle("Delta(1920)", 1920.0, 4, F_Del, FERMION),
	Particle("Delta(1930)", 1950.0, 6, F_Del, FERMION),
	Particle("Delta(1950)", 1930.0, 8, F_Del, FERMION),

	Particle("Delta(2200)", 2200.0, 8, F_Del, FERMION), //!
	Particle("Delta(2420)", 2400.0, 12, F_Del, FERMION),

	//-----------------------------------------------------------------------------
	//Lambda Baryons (pdg starting on pg 83)

	Particle("Lambda", 	  1115.68, 2, F_L, FERMION),
	Particle("Lambda(1405)", 1405.1, 2, F_L, FERMION),
	Particle("Lambda(1520)", 1519.0, 4, F_L, FERMION),
	Particle("Lambda(1600)", 1600.0, 2, F_L, FERMION),
	Particle("Lambda(1670)", 1674.0, 2, F_L, FERMION),
	Particle("Lambda(1690)", 1690.0, 4, F_L, FERMION),
	Particle("Lambda(1800)", 1800.0, 2, F_L, FERMION),
	Particle("Lambda(1810)", 1790.0, 2, F_L, FERMION),
	Particle("Lambda(1820)", 1820.0, 6, F_L, FERMION),
	Particle("Lambda(1830)", 1825.0, 6, F_L, FERMION),
	Particle("Lambda(1890)", 1890.0, 4, F_L, FERMION),

	Particle("Lambda(2100)", 2100.0, 8, F_L, FERMION),
	Particle("Lambda(2110)", 2090.0, 6, F_L, FERMION),
	Particle("Lambda(2350)", 2350.0, 10, F_L,FERMION),

	//-----------------------------------------------------------------------------
	//Sigma Baryons (pdg starting on pg 84)

	Particle("Sigma+", 1189.37, 2, {{1,1,-1},{-1,-1,1}}, FERMION),
	Particle("Sigma0", 1192.64, 2, {{1,0,-1},{-1,0,1}}, FERMION),
	Particle("Sigma-", 1197.45, 2, {{1,-1,-1},{-1,1,1}}, FERMION),
	Particle("Sigma(1385)+", 1382.8, 4, {{1,1,-1},{-1,-1,1}}, FERMION),
	Particle("Sigma(1385)0", 1383.7, 4, {{1,0,-1},{-1,0,1}}, FERMION),
	Particle("Sigma(1385)-", 1387.2, 4, {{1,-1,-1},{-1,1,1}}, FERMION),
	Particle("Sigma(1660)", 1660.0, 2, F_S, FERMION),
	Particle("Sigma(1670)", 1675.0, 4, F_S, FERMION),
	Particle("Sigma(1750)", 1750.0, 2, F_S, FERMION),
	Particle("Sigma(1775)", 1775.0, 6, F_S, FERMION),
	Particle("Sigma(1910)", 1910.0, 4, F_S, FERMION),
	Particle("Sigma(1915)", 1915.0, 6, F_S, FERMION),

	Particle("Sigma(2030)", 2030.0, 8, F_S, FERMION),

	//Excluded from 2022 pdg Summary table, but was in 2012
	// Particle("Sigma(2250)", 2250.0, 2, F_S, FERMION), //J is not listed

	//-----------------------------------------------------------------------------
	// Xi baryons (cascades) (pdg starting on pg 86)

	Particle("Xi0", 1314.86, 2, {{1,0,-2},{-1,0,2}}, FERMION),
	Particle("Xi-", 1321.71, 2, {{1,-1,-2},{-1,1,2}}, FERMION),
	Particle("Xi(1530)0", 1531.8, 4, {{1,0,-2},{-1,0,2}}, FERMION),
	Particle("Xi(1530)-", 1535.0, 4, {{1,-1,-2},{-1,1,2}}, FERMION),
	Particle("Xi(1690)", 1690.0, 2, F_X, FERMION),//J is not listed
	Particle("Xi(1820)", 1823, 4, F_X, FERMION),
	Particle("Xi(1950)", 1950, 2, F_X, FERMION), //J is not listed

	//J is not certain, but it is suggested J >= 5/2
	Particle("Xi(2030)", 2025.0, 6, F_X, FERMION),

	//-----------------------------------------------------------------------------
	//Omega baryons

	Particle("Omega-", 1672.5, 4, F_O, FERMION),

	Particle("Omega(2012)-", 2012.4, 2, F_O, FERMION), //J is not listed //!
	Particle("Omega(2250)-", 2252.0, 2, F_O, FERMION), //J is not listed

	/* Exclude charmed baryons
	//-----------------------------------------------------------------------------
	//Charmed baryons

	Particle("Lambdac+", 2286.5, 2, F_Lc, FERMION),
	Particle("Lambdac(2595)+", 2592.3, 2, F_Lc,FERMION),
	Particle("Lambdac(2625)+", 2628.1, 4, F_Lc,FERMION),
	Particle("Lambdac(2860)+", 2856.1, 4, F_Lc,FERMION),
	Particle("Lambdac(2880)+", 2881.6, 6, F_Lc,FERMION),

	//J is uncertain, using suggested value J=3/2
	Particle("Lambdac(2940)+", 2939.6, 3, F_Lc,FERMION),

	Particle("Sigmac(2455)++", 2454.0, 2, {{1,2,0,1},{-1,-2,0,-1}},FERMION),
	Particle("Sigmac(2455)+", 2452.7, 2, {{1,1,0,1},{-1,-1,0,-1}}, FERMION),
	Particle("Sigmac(2455)0", 2453.8, 2, {{1,0,0,1},{-1,0,0,-1}}, FERMION),
	Particle("Sigmac(2520)++", 2518.4, 4, {{1,2,0,1},{-1,-2,0,-1}},FERMION),
	Particle("Sigmac(2520)+", 2517.4, 4, {{1,1,0,1},{-1,-1,0,-1}}, FERMION,
	Particle("Sigmac(2520)0", 2518.5, 4, {{1,0,0,1},{-1,0,0,-1}}, FERMION),
	Particle("Sigmac(2800)++", 2801.0, 2, {{1,2,0,1},{-1,-2,0,-1}},FERMION),//J is not listed
	Particle("Sigmac(2800)+", 2792.0, 2, {{1,1,0,1},{-1,-1,0,-1}}, FERMION),//J is not listed
	Particle("Sigmac(2800)0", 2806.0, 2, {{1,0,0,1},{-1,0,0,-1}}, FERMION),//J is not listed

	Particle("Xic+", 2467.7, 2, F_Xc1, FERMION),
	Particle("Xic0", 2470.4, 2, F_Xc0, FERMION),
	Particle("Xic'+", 2578.2, 2, F_Xc1, FERMION),
	Particle("Xic'0", 2578.7, 2, F_Xc0, FERMION),

	Particle("Xic(2645)+", 2645.1, 4, F_Xc1, FERMION),
	Particle("Xic(2645)0", 2646.1, 4, F_Xc0, FERMION),
	Particle("Xic(2790)+", 2791.9, 2, F_Xc1, FERMION),
	Particle("Xic(2790)0", 2793.9, 2, F_Xc0, FERMION),
	Particle("Xic(2815)+", 2816.5, 4, F_Xc1, FERMION),
	Particle("Xic(2815)0", 2819.8, 4, F_Xc0, FERMION),
	Particle("Xic(2970)+", 2964.3, 2, F_Xc1, FERMION),
	Particle("Xic(2970)0", 2967.1, 2, F_Xc0, FERMION),

	//J is not listed
	Particle("Xic(3055)", 3055.9, 2,
	 {{1,1,-1,1},{-1,-1,1,-1},{1,0,-1,1},{-1,0,1,-1}}, FERMION), //!

	Particle("Xic(3080)+", 3077.2, 2, F_Xc1, FERMION), //J is not listed
	Particle("Xic(3080)0", 3079.9, 2, F_Xc0, FERMION), //J is not listed

	// Doubly charmed Xi
	Particle("Xicc(3120)++", 3621.6, 2, F_Xcc, FERMION), //J is not listed

	Particle("Omegac0", 2695.2, 2, F_Oc, FERMION),
	Particle("Omegac(2770)0", 2765.9, 4, F_Oc,FERMION),
	Particle("Omegac(3000)0", 3000.4, 2, F_Oc,FERMION), //J is not listed //!
	Particle("Omegac(3050)0", 3050.2, 2, F_Oc,FERMION), //J is not listed //!
	Particle("Omegac(3065)0", 3065.5, 2, F_Oc,FERMION), //J is not listed //!
	Particle("Omegac(3090)0", 3090.1, 2, F_Oc,FERMION), //J is not listed //!
	Particle("Omegac(3120)0", 3119.1, 2, F_Oc,FERMION), //J is not listed //!
	*/
};

#endif