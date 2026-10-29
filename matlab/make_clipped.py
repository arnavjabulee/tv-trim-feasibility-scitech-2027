import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

INPUT_FILE = "aero_pla_127.csv"
LOWER_OUT = "aero_pla_lower.csv"
UPPER_OUT = "aero_pla_upper.csv"

# how aggressive the jump detection is
JUMP_THRESHOLD_MULT = 5.0


def load_xy_csv(fname):
    df = pd.read_csv(fname)
    df = df[["Vt_mps", "Altitude_m"]].dropna().reset_index(drop=True)
    return df


def split_by_jumps(df, thresh_mult=5.0):
    pts = df[["Vt_mps", "Altitude_m"]].to_numpy(dtype=float)

    if len(pts) < 2:
        raise RuntimeError("Not enough points to split.")

    d = np.sqrt(np.sum(np.diff(pts, axis=0) ** 2, axis=1))
    thresh = thresh_mult * np.median(d)
    jump_idx = np.where(d > thresh)[0]

    if len(jump_idx) == 0:
        raise RuntimeError(
            "No large jump found. Try lowering JUMP_THRESHOLD_MULT."
        )

    # use the largest jump if there are multiple
    j = jump_idx[np.argmax(d[jump_idx])]

    lower = df.iloc[: j + 1].copy().reset_index(drop=True)
    upper = df.iloc[j + 1 :].copy().reset_index(drop=True)

    return lower, upper, d, thresh, jump_idx, j


def maybe_reverse_lower(df):
    pts = df[["Vt_mps", "Altitude_m"]].to_numpy()
    # want lower branch to run from low altitude upward
    if pts[0, 1] > pts[-1, 1]:
        df = df.iloc[::-1].reset_index(drop=True)
    return df


def maybe_reverse_upper(df):
    pts = df[["Vt_mps", "Altitude_m"]].to_numpy()
    # want upper branch to run left-to-right in Vt if possible
    if pts[0, 0] > pts[-1, 0]:
        df = df.iloc[::-1].reset_index(drop=True)
    return df


def main():
    df = load_xy_csv(INPUT_FILE)
    lower, upper, d, thresh, jump_idx, split_j = split_by_jumps(
        df, thresh_mult=JUMP_THRESHOLD_MULT
    )

    lower = maybe_reverse_lower(lower)
    upper = maybe_reverse_upper(upper)

    lower.to_csv(LOWER_OUT, index=False)
    upper.to_csv(UPPER_OUT, index=False)

    print(f"Loaded {INPUT_FILE} with {len(df)} points")
    print(f"Median step size: {np.median(d):.6f}")
    print(f"Jump threshold:   {thresh:.6f}")
    print(f"Jump indices:     {jump_idx.tolist()}")
    print(f"Using split at:   {split_j} -> {split_j + 1}")
    print(f"Wrote {LOWER_OUT} ({len(lower)} points)")
    print(f"Wrote {UPPER_OUT} ({len(upper)} points)")

    # debug plots
    plt.figure(figsize=(10, 4))
    plt.plot(d)
    plt.axhline(thresh, linestyle="--")
    plt.xlabel("segment index")
    plt.ylabel("consecutive-point distance")
    plt.title("Jump detection")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("aero_split_jump_debug.png", dpi=160)
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.plot(df["Vt_mps"], df["Altitude_m"], label="raw aero_pla_127")
    plt.xlabel("Vt_mps")
    plt.ylabel("Altitude_m")
    plt.title("Raw curve")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("aero_raw_debug.png", dpi=160)
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.plot(lower["Vt_mps"], lower["Altitude_m"], label="lower")
    plt.plot(upper["Vt_mps"], upper["Altitude_m"], label="upper")
    plt.xlabel("Vt_mps")
    plt.ylabel("Altitude_m")
    plt.title("Split lower / upper")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("aero_split_debug.png", dpi=160)
    plt.close()

    print("Wrote aero_split_jump_debug.png, aero_raw_debug.png, aero_split_debug.png")


if __name__ == "__main__":
    main()