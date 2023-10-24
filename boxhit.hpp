int GetIntersection( float fDst1, float fDst2, CoordFloat P1, CoordPixel P2, CoordFloat& Hit ) {
    if ((fDst1 * fDst2) >= 0.0f) return 0;
    if (fDst1 == fDst2) return 0;
    Hit.x = P1.x + (P2.x - P1.x) * (-fDst1 / (fDst2 - fDst1));
    Hit.y = P1.y + (P2.y - P1.y) * (-fDst1 / (fDst2 - fDst1));
    Hit.z = P1.z + (P2.z - P1.z) * (-fDst1 / (fDst2 - fDst1));
    return 1;
}

int InBox( CoordFloat Hit, CoordFloat B1, CoordFloat B2, const int Axis ) {
    if (Axis == 1 && Hit.z > B1.z && Hit.z < B2.z && Hit.y > B1.y && Hit.y < B2.y) return 1;
    if (Axis == 2 && Hit.z > B1.z && Hit.z < B2.z && Hit.x > B1.x && Hit.x < B2.x) return 1;
    if (Axis == 3 && Hit.x > B1.x && Hit.x < B2.x && Hit.y > B1.y && Hit.y < B2.y) return 1;
    return 0;
}

int CheckLineBox( CoordFloat B1, CoordFloat B2, CoordFloat L1, CoordPixel L2, CoordFloat& Hit1, CoordFloat& Hit2 ) {
    if (L2.x < B1.x && L1.x < B1.x) return false;
    if (L2.x > B2.x && L1.x > B2.x) return false;
    if (L2.y < B1.y && L1.y < B1.y) return false;
    if (L2.y > B2.y && L1.y > B2.y) return false;
    if (L2.z < B1.z && L1.z < B1.z) return false;
    if (L2.z > B2.z && L1.z > B2.z) return false;

    int count = 0;

    CoordFloat TempHit1, TempHit2;

    if (L1.x > B1.x && L1.x < B2.x &&
        L1.y > B1.y && L1.y < B2.y &&
        L1.z > B1.z && L1.z < B2.z) {
        Hit1 = L1;
        count++;
    }

    if (GetIntersection(L1.x - B1.x, L2.x - B1.x, L1, L2, TempHit1)) {
        if (InBox(TempHit1, B1, B2, 1)) {
            if (count == 0) {
                Hit1 = TempHit1;
                count++;
            } else {
                Hit2 = TempHit1;
                return true;
            }
        }
    }

    if (GetIntersection(L1.y - B1.y, L2.y - B1.y, L1, L2, TempHit1)) {
        if (InBox(TempHit1, B1, B2, 2)) {
            if (count == 0) {
                Hit1 = TempHit1;
                count++;
            } else {
                Hit2 = TempHit1;
                return true;
            }
        }
    }

    if (GetIntersection(L1.z - B1.z, L2.z - B1.z, L1, L2, TempHit1)) {
        if (InBox(TempHit1, B1, B2, 3)) {
            if (count == 0) {
                Hit1 = TempHit1;
                count++;
            } else {
                Hit2 = TempHit1;
                return true;
            }
        }
    }

    if (GetIntersection(L1.x - B2.x, L2.x - B2.x, L1, L2, TempHit2)) {
        if (InBox(TempHit2, B1, B2, 1)) {
            if (count == 0) {
                Hit1 = TempHit2;
                count++;
            } else {
                Hit2 = TempHit2;
                return true;
            }
        }
    }

    if (GetIntersection(L1.y - B2.y, L2.y - B2.y, L1, L2, TempHit2)) {
        if (InBox(TempHit2, B1, B2, 2)) {
            if (count == 0) {
                Hit1 = TempHit2;
                count++;
            } else {
                Hit2 = TempHit2;
                return true;
            }
        }
    }

    if (GetIntersection(L1.z - B2.z, L2.z - B2.z, L1, L2, TempHit2)) {
        if (InBox(TempHit2, B1, B2, 3)) {
            if (count == 0) {
                Hit1 = TempHit2;
                count++;
            } else {
                Hit2 = TempHit2;
                return true;
            }
        }
    }

    return (count > 0);
}
