// SignedLog struct and its operations for log-space arithmetic
#[derive(Copy, Clone)]
pub struct SignedLog {
    pub sign: i8, // 1 or -1
    pub log_value: f64,
}

impl SignedLog {
    pub fn log_sum(a: &SignedLog, b: &SignedLog) -> SignedLog {
        if a.sign == 0 { *b }
        else if b.sign == 0 { *a }
        else if a.sign == b.sign {
            let log_value = if a.log_value > b.log_value {
                a.log_value + (1.0 + (b.log_value - a.log_value).exp()).ln()
            } else {
                b.log_value + (1.0 + (a.log_value - b.log_value).exp()).ln()
            };
            SignedLog { sign: a.sign, log_value }
        } else {
            if a.log_value > b.log_value {
                let log_value = a.log_value + (1.0 - (b.log_value - a.log_value).exp()).ln();
                let sign = if log_value.is_infinite() && log_value.is_sign_negative() { 0 } else { a.sign };
                SignedLog { sign, log_value }
            } else {
                let log_value = b.log_value + (1.0 - (a.log_value - b.log_value).exp()).ln();
                let sign = if log_value.is_infinite() && log_value.is_sign_negative() { 0 } else { b.sign };
                SignedLog { sign, log_value }
            }
        }
    }
}