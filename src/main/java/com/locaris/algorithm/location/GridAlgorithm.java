package com.locaris.algorithm.location;

//import java.util.Arrays;
import java.util.ArrayList;

//import org.ejml.data.Complex64F;
//import org.ejml.data.DenseMatrix64F;
//import org.ejml.factory.DecompositionFactory;
//import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.simple.SimpleMatrix;

import com.locaris.algorithm.beans.Position;
import com.locaris.algorithm.tools.CommonFuns;

public class GridAlgorithm implements GridAlgorithmInterface {
	// rdif，一列据差值
	// devicePosition 第一列：X坐标
	// 第二列：Y坐标
	// 第三列：Z坐标
	// 第四列：平方和

	/**********************************************************************************************************
	 ** 函数名字 : pmfTDOA 功能描述 : 基于pmf的TDOA定位算法 作者 : linlin.ma 输 入 : rdif 距差值 :
	 * devicePosition 基站坐标值 : z_hat 标签预置的高度 ： t 连续两次定位间隔 ： position
	 * 上次定位结果：标签坐标及网格信息 输 出 : Position 本次定位结果：标签坐标及网格信息
	 * 
	 *********************************************************************************************************/
	// 网格算法初始化
	@Override
	public Position pmfINITIALIZE(double[][] rdif, double[][] devicePosition, double z_hat) {

		// 位置点
		Position position = new Position();
		try {
			int n;
			n = devicePosition.length;
			SimpleMatrix Pbs = new SimpleMatrix(devicePosition);
			double[] xt;
			double[] yt;

			xt = new double[n];
			yt = new double[n];
			for (int i = 0; i < n; i++) {

				xt[i] = Pbs.get(i, 0);
				yt[i] = Pbs.get(i, 1);
			}

			double x_min = CommonFuns.getMin(xt);
			double x_max = CommonFuns.getMax(xt);
			double y_min = CommonFuns.getMin(yt);
			double y_max = CommonFuns.getMax(yt);

			position.setGridx_min(x_min);
			position.setGridx_max(x_max);
			position.setGridy_min(y_min);
			position.setGridy_max(y_max);

			// 初始化定位结果
			double TagX = (x_min + x_max) / 2;
			double TagY = (y_min + y_max) / 2;
			position.setX(TagX);
			position.setY(TagY);
			return position;
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
			return new Position();
		}
	}

	// 网格算法
	public Position pmfTDOA(double[][] rdif, double[][] devicePosition, double z_hat, double t, Position positionOld) {

		// 位置点
		Position position = new Position();
		// 开始对坐标进行计算
		try {

			SimpleMatrix Rdif = new SimpleMatrix(rdif);
			SimpleMatrix Pbs = new SimpleMatrix(devicePosition);

			double[] tagPos;
			double[] gridBound;

			tagPos = new double[2];
			gridBound = new double[4];

			tagPos[0] = positionOld.getX();
			tagPos[1] = positionOld.getY();
			gridBound[0] = positionOld.getGridx_min();
			gridBound[1] = positionOld.getGridx_max();
			gridBound[2] = positionOld.getGridy_min();
			gridBound[3] = positionOld.getGridy_max();

			double v = 2;
			double noise = 2;
			int n_rd = rdif.length;

			double[] Q = { 2, 0, 0, 2 };
			SimpleMatrix Qm = new SimpleMatrix(2, 2, true, Q);
			SimpleMatrix R = new SimpleMatrix(n_rd, n_rd);
			for (int i = 0; i < n_rd; i++) {
				for (int j = 0; j < n_rd; j++) {
					if (i == j) {
						R.set(i, j, noise + noise);
					} else {
						R.set(i, j, noise);
					}
				}
			}
			SimpleMatrix InvR = R.invert();
			SimpleMatrix InvQ = Qm.invert();
			double Nrm2R = R.conditionP2();
			double Nrm2Qm = Qm.conditionP2();
			double tmp1 = 1 / Math.sqrt(Math.pow(2 * 3.14, n_rd) * Nrm2R);
			double tmp2 = 1 / Math.sqrt(Math.pow(2 * 3.14, 2) * Nrm2Qm);

			// 时间更新，网格边界扩张
			double[] tagPostmp;
			double[] gridBoundtmp;
			tagPostmp = new double[2];
			gridBoundtmp = new double[4];
			gridBoundtmp[0] = gridBound[0] - v * t;
			gridBoundtmp[1] = gridBound[1] + v * t;
			gridBoundtmp[2] = gridBound[2] - v * t;
			gridBoundtmp[3] = gridBound[3] + v * t;

			// 量测更新
			// 第一步粗化网格
			double res_first = 0.2;
			double res_second = 0.1;
			double thr = 0.01;
			int Num, xnum, ynum;
			double Lx, Ly;
			double[] ref_grid;
			ref_grid = new double[2];

			ref_grid[0] = gridBoundtmp[0];
			ref_grid[1] = gridBoundtmp[2];
			Lx = gridBoundtmp[1] - gridBoundtmp[0];
			Ly = gridBoundtmp[3] - gridBoundtmp[2];
			xnum = (int) (Lx / res_first + 1 + 0.5);
			ynum = (int) (Ly / res_first + 1 + 0.5);
			Num = xnum * ynum;
			// 网格生成
			SimpleMatrix gridx = new SimpleMatrix(Num, 1);
			SimpleMatrix gridy = new SimpleMatrix(Num, 1);
			for (int i = 1; i <= ynum; i++) {
				for (int j = 1; j <= xnum; j++) {
					gridx.set((i - 1) * xnum + j - 1, 0, (j - 1) * res_first + ref_grid[0]);
					gridy.set((i - 1) * xnum + j - 1, 0, (i - 1) * res_first + ref_grid[1]);
				}
			}
			System.out.println("num:"+Num);
			SimpleMatrix diff_distance = new SimpleMatrix(n_rd, Num);// 网格点到相应基站的距差
			SimpleMatrix ref_distance = new SimpleMatrix(Num, 1); // 网格点到参考基站的距离
			SimpleMatrix res_err = new SimpleMatrix(Num, 1); // 网格点距差和量测距差的误差系数
			SimpleMatrix piror_pdf = new SimpleMatrix(Num, 1); // 网格点的先验概率
			SimpleMatrix post_pdf = new SimpleMatrix(Num, 1); // 网格点的后验概率
			double rx, ry, rz;
			for (int i = 0; i < Num; i++) {
				piror_pdf.set(i, 0, 1.0 / Num);
				rx = (Pbs.get(0, 0) - gridx.get(i, 0)) * (Pbs.get(0, 0) - gridx.get(i, 0));
				ry = (Pbs.get(0, 1) - gridy.get(i, 0)) * (Pbs.get(0, 1) - gridy.get(i, 0));
				rz = (Pbs.get(0, 2) - z_hat) * (Pbs.get(0, 2) - z_hat);
				ref_distance.set(i, 0, Math.sqrt(rx + ry + rz)); // 计算网格点到参考基站的距离
			}
			SimpleMatrix dx = new SimpleMatrix(n_rd, 1);
			SimpleMatrix dy = new SimpleMatrix(n_rd, 1);
			SimpleMatrix dz = new SimpleMatrix(n_rd, 1);
			double[] dis_tmp, res;
			dis_tmp = new double[n_rd];
			res = new double[n_rd];
			SimpleMatrix res1 = new SimpleMatrix(1, n_rd);
			SimpleMatrix res2 = new SimpleMatrix(1, 1);
			double sum = 0;
			for (int i = 0; i < Num; i++) {
				for (int k = 0; k < n_rd; k++) {
					dx.set(k, 0, Pbs.get(k + 1, 0) - gridx.get(i, 0));
					dy.set(k, 0, Pbs.get(k + 1, 1) - gridy.get(i, 0));
					dz.set(k, 0, Pbs.get(k + 1, 2) - z_hat);
				}
				SimpleMatrix dxm = dx.elementMult(dx);
				SimpleMatrix dym = dy.elementMult(dy);
				SimpleMatrix dzm = dz.elementMult(dz);
				for (int j = 0; j < n_rd; j++) {
					dis_tmp[j] = Math.sqrt(dxm.get(j, 0) + dym.get(j, 0) + dzm.get(j, 0)) - ref_distance.get(i, 0);
				}
				diff_distance.setColumn(i, 0, dis_tmp);
				for (int j = 0; j < n_rd; j++) {
					res[j] = diff_distance.get(j, i) - Rdif.get(j, 0);
				}
				SimpleMatrix res_tmp = new SimpleMatrix(n_rd, 1, true, res);
				res1 = res_tmp.transpose().mult(InvR);
				res2 = res1.mult(res_tmp);
				res_err.set(i, 0, res2.get(0, 0));
				post_pdf.set(i, 0, tmp1 * Math.exp(-1.0 / 2.0 * res_err.get(i, 0)) * piror_pdf.get(i, 0));
				sum = sum + post_pdf.get(i, 0);
			}
			double[] reserrtmp;
			reserrtmp = new double[Num];
			double reserr_min;
			for (int i = 0; i < Num; i++) {
				reserrtmp[i] = res_err.get(i, 0);
			}
			reserr_min = CommonFuns.getMin(reserrtmp);
			if (reserr_min > 1) {
				gridBoundtmp[0] = gridBoundtmp[0] - 1;
				gridBoundtmp[1] = gridBoundtmp[1] + 1;
				gridBoundtmp[2] = gridBoundtmp[2] - 1;
				gridBoundtmp[3] = gridBoundtmp[3] + 1;
				ref_grid[0] = gridBoundtmp[0];
				ref_grid[1] = gridBoundtmp[2];
				Lx = gridBoundtmp[1] - gridBoundtmp[0];
				Ly = gridBoundtmp[3] - gridBoundtmp[2];
				xnum = (int) (Lx / res_first + 1 + 0.5);
				ynum = (int) (Ly / res_first + 1 + 0.5);
				Num = xnum * ynum;
				// 网格扩张后，网格点数增大，需要对矩阵维度进行改变，需要重新定义大小
				gridx.reshape(Num, 1);
				gridy.reshape(Num, 1);
				for (int i = 1; i <= ynum; i++) {
					for (int j = 1; j <= xnum; j++) {
						gridx.set((i - 1) * xnum + j - 1, 0, (j - 1) * res_first + ref_grid[0]);
						gridy.set((i - 1) * xnum + j - 1, 0, (i - 1) * res_first + ref_grid[1]);
					}
				}
				diff_distance.reshape(n_rd, Num);
				ref_distance.reshape(Num, 1);
				res_err.reshape(Num, 1);
				piror_pdf.reshape(Num, 1);
				post_pdf.reshape(Num, 1);
				for (int i = 0; i < Num; i++) {
					piror_pdf.set(i, 0, 1.0 / Num);
					rx = (Pbs.get(0, 0) - gridx.get(i, 0)) * (Pbs.get(0, 0) - gridx.get(i, 0));
					ry = (Pbs.get(0, 1) - gridy.get(i, 0)) * (Pbs.get(0, 1) - gridy.get(i, 0));
					rz = (Pbs.get(0, 2) - z_hat) * (Pbs.get(0, 2) - z_hat);
					ref_distance.set(i, 0, Math.sqrt(rx + ry + rz));
				}
				sum = 0;
				dx.reshape(n_rd, 1);
				dy.reshape(n_rd, 1);
				dz.reshape(n_rd, 1);
				res1.reshape(1, n_rd);
				res2.reshape(1, 1);
				for (int i = 0; i < Num; i++) {
					for (int k = 0; k < n_rd; k++) {
						dx.set(k, 0, Pbs.get(k + 1, 0) - gridx.get(i, 0));
						dy.set(k, 0, Pbs.get(k + 1, 1) - gridy.get(i, 0));
						dz.set(k, 0, Pbs.get(k + 1, 2) - z_hat);
					}
					SimpleMatrix dxm = dx.elementMult(dx);
					SimpleMatrix dym = dy.elementMult(dy);
					SimpleMatrix dzm = dz.elementMult(dz);
					for (int j = 0; j < n_rd; j++) {
						dis_tmp[j] = Math.sqrt(dxm.get(j, 0) + dym.get(j, 0) + dzm.get(j, 0)) - ref_distance.get(i, 0);
					}
					diff_distance.setColumn(i, 0, dis_tmp);
					for (int j = 0; j < n_rd; j++) {
						res[j] = diff_distance.get(j, i) - Rdif.get(j, 0);
					}
					SimpleMatrix res_tmp = new SimpleMatrix(n_rd, 1, true, res);
					res1 = res_tmp.transpose().mult(InvR);
					res2 = res1.mult(res_tmp);
					res_err.set(i, 0, res2.get(0, 0));
					post_pdf.set(i, 0, tmp1 * Math.exp(-1 / 2 * res_err.get(i, 0)) * piror_pdf.get(i, 0));
					sum = sum + post_pdf.get(i, 0);
				}
				double themp1 = thr / Num / (res_first * res_first);
				if (sum == 0) {
					tagPostmp[0] = tagPos[0];
					tagPostmp[1] = tagPos[1];
					gridBoundtmp[0] = gridBoundtmp[0] - 2;
					gridBoundtmp[1] = gridBoundtmp[1] + 2;
					gridBoundtmp[2] = gridBoundtmp[2] - 2;
					gridBoundtmp[3] = gridBoundtmp[3] + 2;
				} else {
					for (int i = 0; i < Num; i++) {
						post_pdf.set(i, 0, post_pdf.get(i, 0) / sum);
						tagPostmp[0] = tagPostmp[0] + gridx.get(i, 0) * post_pdf.get(i, 0);
						tagPostmp[1] = tagPostmp[1] + gridy.get(i, 0) * post_pdf.get(i, 0);
					}
					ArrayList<Double> gridxlist = new ArrayList<Double>();
					ArrayList<Double> gridylist = new ArrayList<Double>();
					for (int i = 0; i < Num; i++) {
						if (post_pdf.get(i, 0) >= themp1) {
							gridxlist.add(gridx.get(i, 0));
							gridylist.add(gridy.get(i, 0));
						}
					}
					if (gridxlist.size() <= 4) {
						gridBoundtmp[0] = gridBoundtmp[0] - 1;
						gridBoundtmp[1] = gridBoundtmp[1] + 1;
						gridBoundtmp[2] = gridBoundtmp[2] - 1;
						gridBoundtmp[3] = gridBoundtmp[3] + 1;
					} else {
						Double[] gridxrem = new Double[gridxlist.size()];
						Double[] gridyrem = new Double[gridylist.size()];
						gridxrem = gridxlist.toArray(gridxrem);
						gridyrem = gridylist.toArray(gridyrem);
						gridBoundtmp[0] = CommonFuns.getMin(gridxrem);
						gridBoundtmp[1] = CommonFuns.getMax(gridxrem);
						gridBoundtmp[2] = CommonFuns.getMin(gridyrem);
						gridBoundtmp[3] = CommonFuns.getMax(gridyrem);
					}
				}
			} else {
				double themp1 = thr / Num / (res_first * res_first);
				if (sum == 0) {
					tagPostmp[0] = tagPos[0];
					tagPostmp[1] = tagPos[1];
					gridBoundtmp[0] = gridBoundtmp[0] - 2;
					gridBoundtmp[1] = gridBoundtmp[1] + 2;
					gridBoundtmp[2] = gridBoundtmp[2] - 2;
					gridBoundtmp[3] = gridBoundtmp[3] + 2;
				} else {
					for (int i = 0; i < Num; i++) {
						post_pdf.set(i, 0, post_pdf.get(i, 0) / sum);
						tagPostmp[0] = tagPostmp[0] + gridx.get(i, 0) * post_pdf.get(i, 0);
						tagPostmp[1] = tagPostmp[1] + gridy.get(i, 0) * post_pdf.get(i, 0);
					}
					ArrayList<Double> gridxlist = new ArrayList<Double>();
					ArrayList<Double> gridylist = new ArrayList<Double>();
					for (int i = 0; i < Num; i++) {
						if (post_pdf.get(i, 0) >= themp1) {
							gridxlist.add(gridx.get(i, 0));
							gridylist.add(gridy.get(i, 0));
						}
					}
					if (gridxlist.size() <= 4) {
						gridBoundtmp[0] = gridBoundtmp[0] - 1;
						gridBoundtmp[1] = gridBoundtmp[1] + 1;
						gridBoundtmp[2] = gridBoundtmp[2] - 1;
						gridBoundtmp[3] = gridBoundtmp[3] + 1;
					} else {
						Double[] gridxrem = new Double[gridxlist.size()];
						Double[] gridyrem = new Double[gridylist.size()];
						gridxrem = gridxlist.toArray(gridxrem);
						gridyrem = gridylist.toArray(gridyrem);
						gridBoundtmp[0] = CommonFuns.getMin(gridxrem);
						gridBoundtmp[1] = CommonFuns.getMax(gridxrem);
						gridBoundtmp[2] = CommonFuns.getMin(gridyrem);
						gridBoundtmp[3] = CommonFuns.getMax(gridyrem);
					}
					// 第二步细化网格
					ref_grid[0] = gridBoundtmp[0];
					ref_grid[1] = gridBoundtmp[2];
					Lx = gridBoundtmp[1] - gridBoundtmp[0];
					Ly = gridBoundtmp[3] - gridBoundtmp[2];
					xnum = (int) (Lx / res_second + 1 + 0.5);
					ynum = (int) (Ly / res_second + 1 + 0.5);
					Num = xnum * ynum;
					// 网格生成
					SimpleMatrix grid = new SimpleMatrix(Num, 2);
					for (int i = 1; i <= ynum; i++) {
						for (int j = 1; j <= xnum; j++) {
							grid.set((i - 1) * xnum + j - 1, 0, (j - 1) * res_second + ref_grid[0]);
							grid.set((i - 1) * xnum + j - 1, 1, (i - 1) * res_second + ref_grid[1]);
						}
					}
					SimpleMatrix xyPos = new SimpleMatrix(1, 2, true, tagPostmp);
					SimpleMatrix piror_pdf2 = new SimpleMatrix(Num, 1);
					SimpleMatrix p1 = new SimpleMatrix(1, 2);
					SimpleMatrix p2 = new SimpleMatrix(1, 1);
					SimpleMatrix xytmp = new SimpleMatrix(1, 2);
					SimpleMatrix xytmpT = new SimpleMatrix(2, 1);
					for (int i = 0; i < Num; i++) {
						xytmp.set(0, 0, grid.get(i, 0) - xyPos.get(0, 0));
						xytmp.set(0, 1, grid.get(i, 1) - xyPos.get(0, 1));
						xytmpT.reshape(2, 1);
						p1.reshape(1, 2);
						p2.reshape(1, 1);
						xytmpT = xytmp.transpose();
						p1 = xytmp.mult(InvQ);
						p2 = p1.mult(xytmpT);
						piror_pdf2.set(i, 0, tmp2 * Math.exp(-1.0 / 2.0 * p2.get(0, 0)));
					}
					SimpleMatrix diff_distance2 = new SimpleMatrix(n_rd, Num);
					SimpleMatrix ref_distance2 = new SimpleMatrix(Num, 1);
					SimpleMatrix res_err2 = new SimpleMatrix(Num, 1);
					SimpleMatrix post_pdf2 = new SimpleMatrix(Num, 1);
					for (int i = 0; i < Num; i++) {
						rx = (Pbs.get(0, 0) - grid.get(i, 0)) * (Pbs.get(0, 0) - grid.get(i, 0));
						ry = (Pbs.get(0, 1) - grid.get(i, 1)) * (Pbs.get(0, 1) - grid.get(i, 1));
						rz = (Pbs.get(0, 2) - z_hat) * (Pbs.get(0, 2) - z_hat);
						ref_distance2.set(i, 0, Math.sqrt(rx + ry + rz));
					}
					sum = 0;
					for (int i = 0; i < Num; i++) {
						for (int k = 0; k < n_rd; k++) {
							dx.set(k, 0, Pbs.get(k + 1, 0) - grid.get(i, 0));
							dy.set(k, 0, Pbs.get(k + 1, 1) - grid.get(i, 1));
							dz.set(k, 0, Pbs.get(k + 1, 2) - z_hat);
						}
						SimpleMatrix dxm = dx.elementMult(dx);
						SimpleMatrix dym = dy.elementMult(dy);
						SimpleMatrix dzm = dz.elementMult(dz);
						for (int j = 0; j < n_rd; j++) {
							dis_tmp[j] = Math.sqrt(dxm.get(j, 0) + dym.get(j, 0) + dzm.get(j, 0))
									- ref_distance2.get(i, 0);
						}
						diff_distance2.setColumn(i, 0, dis_tmp);
						for (int j = 0; j < n_rd; j++) {
							res[j] = diff_distance2.get(j, i) - Rdif.get(j, 0);
						}
						SimpleMatrix res_tmp = new SimpleMatrix(n_rd, 1, true, res);
						res1 = res_tmp.transpose().mult(InvR);
						res2 = res1.mult(res_tmp);
						res_err2.set(i, 0, res2.get(0, 0));
						post_pdf2.set(i, 0, tmp1 * Math.exp(-1 / 2 * res_err2.get(i, 0)) * piror_pdf2.get(i, 0));
						sum = sum + post_pdf2.get(i, 0);
					}
					double themp2 = thr / Num / (res_second * res_second);
					if (sum == 0) {
						tagPostmp[0] = tagPostmp[0];
						tagPostmp[1] = tagPostmp[1];
						gridBoundtmp[0] = gridBoundtmp[0];
						gridBoundtmp[1] = gridBoundtmp[1];
						gridBoundtmp[2] = gridBoundtmp[2];
						gridBoundtmp[3] = gridBoundtmp[3];
					} else {
						tagPostmp[0] = 0;
						tagPostmp[1] = 0;
						for (int i = 0; i < Num; i++) {
							post_pdf2.set(i, 0, post_pdf2.get(i, 0) / sum);
							tagPostmp[0] = tagPostmp[0] + grid.get(i, 0) * post_pdf2.get(i, 0);
							tagPostmp[1] = tagPostmp[1] + grid.get(i, 1) * post_pdf2.get(i, 0);
						}
						ArrayList<Double> gridxlist2 = new ArrayList<Double>();
						ArrayList<Double> gridylist2 = new ArrayList<Double>();
						for (int i = 0; i < Num; i++) {
							if (post_pdf2.get(i, 0) >= themp2) {
								gridxlist2.add(grid.get(i, 0));
								gridylist2.add(grid.get(i, 1));
							}
						}
						if (gridxlist2.size() <= 4) {
							gridBoundtmp[0] = gridBoundtmp[0] - 1;
							gridBoundtmp[1] = gridBoundtmp[1] + 1;
							gridBoundtmp[2] = gridBoundtmp[2] - 1;
							gridBoundtmp[3] = gridBoundtmp[3] + 1;
						} else {
							Double[] gridxrem2 = new Double[gridxlist2.size()];
							Double[] gridyrem2 = new Double[gridylist2.size()];
							gridxrem2 = gridxlist2.toArray(gridxrem2);
							gridyrem2 = gridylist2.toArray(gridyrem2);
							gridBoundtmp[0] = CommonFuns.getMin(gridxrem2);
							gridBoundtmp[1] = CommonFuns.getMax(gridxrem2);
							gridBoundtmp[2] = CommonFuns.getMin(gridyrem2);
							gridBoundtmp[3] = CommonFuns.getMax(gridyrem2);
						}
					}
				}
			}
			position.setX(tagPostmp[0]);
			position.setY(tagPostmp[1]);
			position.setGridx_min(gridBoundtmp[0]);
			position.setGridx_max(gridBoundtmp[1]);
			position.setGridy_min(gridBoundtmp[2]);
			position.setGridy_max(gridBoundtmp[3]);
			return position;
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
			return new Position();
		}
	}

}
