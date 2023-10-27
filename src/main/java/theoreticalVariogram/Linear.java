/*
 * GNU GPL v3 License
 *
 * Copyright 2016 Marialaura Bancheri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package theoreticalVariogram;

public class Linear implements Model {

	double dist;
	double sill;
	double range;
	double nug;
	boolean isOk = false;

	public Linear(double dist, double sill, double range, double nug) {
		this.dist = dist;
		this.sill = sill;
		this.range = range;
		this.nug = nug;
		this.isOk = nug >= 0 && sill >= 0 && range >= 0;

	}

	@Override
	public double computeSemivariance() {

		double result = 0;

		if (dist > 0.0 & dist <= range) {
			result = nug + sill * (dist / range);
		}
		if (dist > range) {
			result = sill + nug;
		}

		return result;
	}

	@Override
	public double[] computeGradient() {
		double[] gradient = new double[] { Double.NaN, Double.NaN, Double.NaN };
		if (dist > 0.0 & dist <= range) {
			gradient = new double[] { 1.0, dist / range, -(sill * dist) / Math.pow(range, 2.0) };
		} else if (dist > range) {
			gradient = new double[] { 1.0, 1.0, 0.0 };
		}
		return gradient;
	}

}
