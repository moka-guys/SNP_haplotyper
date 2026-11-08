CREATE TABLE IF NOT EXISTS app_users (
    id                   SERIAL PRIMARY KEY,
    username             VARCHAR(50) UNIQUE NOT NULL,
    password_hash        VARCHAR(255) NOT NULL,
    is_active            BOOLEAN NOT NULL DEFAULT TRUE,
    created_at           TIMESTAMPTZ DEFAULT NOW()
);